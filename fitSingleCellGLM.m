% ff = rdir('D:\Neuroscience\ErasmusAttn\glmExpt_*_25ms.mat');
% for ii = 1:numel(ff); fitSingleCellGLM(ff(ii).name); end
function fitSingleCellGLM(dataFile, baselineTrials, lazy)

  %% Default arguments
  if ~exist('baselineTrials', 'var')
    baselineTrials    = 5;
  end
  if ~exist('lazy', 'var') || isempty(lazy)
    lazy              = true;
  end
  
  %% Setup optimization package
  optimPackage        = dir(fullfile(fileparts(mfilename('fullpath')), 'spams*', 'start_*.m'));
  origDir             = cd(optimPackage.folder);
  run(optimPackage.name);
  cd(origDir);

  %% Analysis configuration
  cfg                 = struct();
  
  % Cross validation and permutation tests
  cfg.nCVFolds        = 5;
  cfg.nLambdas        = 10;
  cfg.minLambda       = 1e-4;
  cfg.maxRelLikeli    = 0.1;            % maximum relative likelihood of models to continue specialization
  
  % Trial selection criteria
  cfg.maxCueStart     = 1000;           % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
%   cfg.selectPastCond  = 0:4;            % keep only trials with these past_condition
  cfg.selectPastCond  = 1:2;            % keep only trials with these past_condition
  cfg.minNumTrials    = cfg.nCVFolds;   % minimum number of trials for fitting models
  
  % Maximum duration of neural responses per behavioral event
  cfg.eventDuration   = 1500;           % in ms
%   cfg.eventDuration   = 1000;           % in ms
  cfg.eventNBasis     = 5;              % number of basis functions to use per event-response
%   cfg.eventDuration   = 2000;           % in ms
%   cfg.eventNBasis     = 5;              % number of basis functions to use per event-response
  
  % Configuration for SPAMS
  cfg.fitOptions      = struct('loss', 'poisson-log', 'regul', 'group-lasso-l2', 'intercept', true, 'max_it', 1e5, 'verbose', false, 'L0', 1e-5, 'ista', true, 'linesearch_mode', 2);

  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
%   name                = sprintf('glmFit_%s_min%dtrials_%.0flikeli', regexprep(name, '^[^_]+_', ''), cfg.minNumTrials, 100*cfg.maxRelLikeli);
  name                = sprintf('glmFit_%s_min%dtrials_%.0fdur', regexprep(name, '^[^_]+_', ''), cfg.minNumTrials, cfg.eventDuration);
  if ~isempty(baselineTrials)
    name              = sprintf('%s_baseline%d', name, baselineTrials);
  end
  outputFile          = fullfile(path, [name ext]);
  
  %% Specify experimental variables that are included in the model
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
  cfg.behavConditions = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  
  % Overly many basis functions can lead to overfitting
  cfg.behavEvents(~cellfun(@isempty, regexp(cfg.behavEvents, '(^cue_|_stop$)', 'once')))  = [];
  
  % Not sure how to deal with these yet
  cfg.behavConditions = setdiff(cfg.behavConditions, {'trial_nr', 'reward_duration', 'saccade_amplitude'});
  
  % Too many categories to search over
%   cfg.behavConditions( ~cellfun(@isempty, regexp(cfg.behavConditions, '^past_(?!condition)', 'once')) ) = [];

  %% Use same trial categories as behavior model
  cfg.behavConditions(cellfun(@isempty, regexp(cfg.behavConditions, 'gap|saccade', 'once')))  = [];
  
  %% Support resuming of previously saved fits
  CategorizedModel.nextID(0);
  if lazy && exist(outputFile, 'file')
    newCfg            = cfg;
    load(outputFile);
    if ~isequal(newCfg, keepfield(cfg,fieldnames(newCfg)))
      error('fitSingleCellGLM:cfg', 'Inconsistent analysis configuration found in old model file %s', outputFile);
    end
  else
    unspecializedModel= cell(size(data.cellData));
    categoryModel     = cell(size(data.cellData));
    hierarchicalModel = cell(size(data.cellData));
  end
  
  
  %% Configure models per cell
  for iCell = 1:numel(data.cellData)
    if ~isempty(unspecializedModel{iCell})
      continue
    end
    
    %% No reward indicator for error trials
    for iTrial = find([data.cellData(iCell).trials.reward_duration] == 0)
      data.cellData(iCell).trials(iTrial).reward_start  = nan;
    end
    
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                      ...
                      & ~any(diff(accumfun(1, @(x) [data.cellData(iCell).trials.(x)], cfg.behavEvents), 1, 1) <= data.experiment.binSize, 1)  ...
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)    ...
                      & ismember([data.cellData(iCell).trials.past_condition], cfg.selectPastCond)      ...
                      & arrayfun(@(x) ~isempty(x.SS), data.cellData(iCell).trials(:)')                  ... not sure why some trials can have zero spikes
                      & arrayfun(@(x) ~isnan(x.saccade_onset), data.cellData(iCell).trials(:)')         ... must have timing info for all events 
                      ;
%     figure; imagesc(diff(accumfun(1, @(x) [data.cellData(iCell).trials.(x)], cfg.behavEvents), 1, 1) > data.experiment.binSize)
                    
    experiment        = data.experiment;
    experiment.id     = [data.cellData(iCell).monkey '_' data.cellData(iCell).cell_id];
    experiment.trial  = data.cellData(iCell).trials(selTrials);
    experiment.selectedTrials   = selTrials;
    
    if numel(experiment.trial) < cfg.minNumTrials
      unspecializedModel{iCell} = experiment;
      continue;
    end
    
%     assert(all( arrayfun(@(x) all(isfinite(cellfun(@(y) x.(y), cfg.behavEvents))), experiment.trial) ));
    
    %% Divide trials into categories according to various stimulus and outcome variables
    trialConditions   = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);
    [trialCategory, ~, categIndex]    = unique(trialConditions, 'rows');
    assert(~any(isnan(trialCategory(:))));
    
    % Construct a joint code for unique data categories
    categDigits       = 1 + max(floor(log10(max(1,trialCategory))), [], 1);
    categCode         = sum(trialCategory .* [10.^cumsum(categDigits(2:end), 'reverse'), 1], 2);
    trialCode         = num2cell(categCode(categIndex));
    [experiment.trial.behavior_code]  = trialCode{:};
    
    %% Specify the form of neural response timecourses to various experimental variables
    modelSpecs        = buildGLM.initDesignSpec(experiment);

    % Responses to behavioral events that vary smoothly with time after the event
    basisFcn          = basisFactory.makeSmoothTemporalBasis('progressive cosine x1.5 off0', cfg.eventDuration, cfg.eventNBasis, experiment.binfun);
    basisFcn.B        = basisFcn.B ./ max(basisFcn.B,[],1);
%     figure; plot(basisFcn.tr * experiment.binSize,basisFcn.B)
    for iVar = 1:numel(cfg.behavEvents)
      modelSpecs      = buildGLM.addCovariateTiming(modelSpecs, cfg.behavEvents{iVar}, cfg.behavEvents{iVar}, experiment.desc.(cfg.behavEvents{iVar}), basisFcn);
    end

    %% Specify design matrix using the above parameters
    [design,trialBin,eventEpochs] = buildGLM.compileSparseDesignMatrix(modelSpecs, 1:numel(experiment.trial));
    design.dspec.expt.numTimeBins = size(design.X,1);
    design.dspec.expt.trialBins   = trialBin;
    design.dspec.expt.eventEpochs = eventEpochs;

    %% Optional scaling of regressors for each trial by the moving average firing rate 
    if ~isempty(baselineTrials)
      meanRate        = arrayfun(@(x) numel(x.SS)/x.duration, experiment.trial);    % rate in 1/ms
%       longfigure; hold on; plot(meanRate); plot(movingAverage(meanRate, baselineTrials, 1))
      meanRate        = movingAverage(meanRate, baselineTrials, 1);
      meanRate        = log(meanRate * experiment.binSize);                         % convert rate to the time bins used
      design.X(:,end+1) = meanRate(design.trialX);
      design.biasCol  = size(design.X,2);
    end

    %% Store model specifications and data for this cell
    firingRate        = full(buildGLM.getBinnedSpikeTrain(experiment, 'SS', design.trialIndices));
    unspecializedModel{iCell}     = CategorizedModel(firingRate, design);
  end
  
  %% Fit model with no specialization of regressors by behavior categories
  fprintf('Fitting unspecialized models for %d/%d cells...\n', sum(cellfun(@(x) ~isstruct(x) && isempty(x.prediction), unspecializedModel)), numel(data.cellData));
  parfor iCell = 1:numel(unspecializedModel)
    if ~isstruct(unspecializedModel{iCell}) && isempty(unspecializedModel{iCell}.prediction)
      unspecializedModel{iCell} = unspecializedModel{iCell}.regress(cfg.nCVFolds, cfg.fitOptions, false, cfg.nLambdas, cfg.minLambda);
    end
  end
  
  %% Progressive save
  save(outputFile, 'unspecializedModel', 'categoryModel', 'hierarchicalModel', 'cfg');
  
  %% Fit categorized and hierarchical models per cell
  fprintf('Fitting specialized models:\n');
  for iCell = 1:numel(unspecializedModel)
    fullModel         = unspecializedModel{iCell};
    if isstruct(fullModel) || ~isempty(hierarchicalModel{iCell})
      continue
    end
    
    tStart            = tic;
    fprintf('  [%3d]  %-20s', iCell, fullModel.design.dspec.expt.id);
    drawnow;
      
    %% Unused : Fit model separately for data per unique category
    %{
    if isempty(categoryModel{iCell})
      categModel      = fullModel.specialize('behavior_code', false, cfg.minNumTrials);
      categModel      = categModel{1};
      categoryModel{iCell}    = categModel.regress(cfg.nCVFolds, cfg.fitOptions, true, cfg.nLambdas, cfg.minLambda);
    end
    %}
  
    %% Hierarchical model to find a parsimonious set of trial-category-based specializations for this neuron's response
    if isempty(hierarchicalModel{iCell})
      hierarchicalModel{iCell}= buildHierarchicalModel(fullModel, cfg.behavConditions, cfg);
    end

    %% Progressive save
    save(outputFile, 'unspecializedModel', 'categoryModel', 'hierarchicalModel', 'cfg');
    fprintf(' ... %gs\n', toc(tStart));
  end
  
  
  fprintf(' -->  %s\n', outputFile);
  
end

function [model, altModels] = buildHierarchicalModel(basicModel, conditionLabels, cfg)
  
  model                     = basicModel;
  altModels                 = {};
  
  %% Continue specializing models until there is no improvement in goodness of fit
  while ~isempty(conditionLabels)
    %% Generate further specialized models per candidate condition 
    candModel               = {};
    for iCond = 1:numel(conditionLabels)
      %% Construct models for all possible partitions of values for a given condition label
      tryModels             = model.specialize(conditionLabels{iCond}, true, cfg.minNumTrials);
%       tryModels( cellfun(@(x) any(x.numberOfTrials() < cfg.minNumTrials), tryModels) )  = [];
      candModel             = [candModel; tryModels];
    end
    if isempty(candModel)
      break;
    end
    
    %% Parallel optimization for all candidate models; N.B. handle objects are modified in place
    allModels               = [candModel{:}];
    allModels.regress(cfg.nCVFolds, cfg.fitOptions, true, cfg.nLambdas, cfg.minLambda);
    
    %% Select model with the best relative likelihood w.r.t. the less specialized parent model
    relLikeli               = cellfun(@(x) x.relativeLikelihood(), candModel);
    if any(isnan(relLikeli))
      warning('buildHierarchicalModel:relLikeli', 'Invalid log-likelihood encountered for model(s) %s', mat2str(find(~isfinite(relLikeli))));
      veto                  = ~isfinite(relLikeli);
      relLikeli(veto)       = [];
      candModel(veto)       = [];
    end
    
    [sortedRelLikeli,iOrder]= sort(relLikeli, 'ascend');
    altModels{end+1}        = candModel(iOrder(2:end));
%     longfigure;   hold on;  [predict,target]=model.jointPrediction; plot(predict);  plot(accumfun(2, @(x) x.jointPrediction, candModel(iOrder([1 end]))));    plot(target,'k');   
    
    %% Check if we should continue specializing models
    if sortedRelLikeli(1) > cfg.maxRelLikeli
      break;
    end
    
    model                   = candModel{iOrder(1)};
    assert(all(strcmp({model(2:end).categoryLabel}, model(1).categoryLabel)));
    conditionLabels( strcmp(conditionLabels,model(1).categoryLabel) ) = [];
    
    %% Diagnostics
    %{
    iCand                   = iOrder(1);
    categValues             = arrayfun(@(x) mat2str(x.categoryValues'), candModel{iCand}, 'UniformOutput', false);
    categLegend             = [strrep(candModel{iCand}(1).categoryLabel,'_',' '), ' = ', strjoin(coloredText(categValues, lines(numel(categValues))), ' / ')];
    [pan,shape]             = makePanels(numel(fieldnames(model.regressors)), '', '', 'maxsubplotsize', 200);
    iPlot                   = 0;
    
    for what = fieldnames(model.regressors)'
      if ~isstruct(model.regressors.(what{:}))
        continue
      end
      
      iPlot                 = iPlot + 1;
      axs                   = selectPanel(pan, iPlot, shape);
      hold(axs, 'on');
      title(axs, categLegend, 'fontweight', 'normal');
      xlabel(axs, sprintf('time from event (%s)', model.design.dspec.expt.unitOfTime));
      ylabel(axs, [strrep(what{:}, '_', ' ') ' response']);
      
      for iModel = 1:numel(candModel{iCand})
        plot(axs, candModel{iCand}(iModel).regressors.(what{:}).time, candModel{iCand}(iModel).regressors.(what{:}).response);
      end
      plot(axs, model.regressors.(what{:}).time, model.regressors.(what{:}).response, 'k');
    end
    %}
  end
  
end
