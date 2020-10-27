function fitSingleCellGLM(dataFile, postfix, lazy)

  %% Default arguments
  if ~exist('postfix', 'var') || isempty(postfix)
    postfix           = '';
  end
  if ~exist('lazy', 'var') || isempty(lazy)
    lazy              = true;
  end

  %% Analysis configuration
  cfg                 = struct();
  
  % Trial selection criteria
  cfg.maxCueStart     = 300;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
  cfg.minNumTrials    = 10;             % minimum number of trials for fitting models
  
  % Maximum duration of neural responses per behavioral event
  cfg.eventDuration   = 1500;           % in ms
  cfg.eventNBasis     = 4;              % number of basis functions to use per event-response
  
  % Cross validation and permutation tests
  cfg.nShuffles       = 100;
  cfg.nBootstrap      = 0;
  cfg.nCVFolds        = 3;
  cfg.nLambdas        = 5;
  cfg.maxRelLikeli    = 0.05;           % maximum relative likelihood of models to continue specialization
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = ['glmFit_' regexprep(name, '^[^_]+_', '')];
  outputFile          = fullfile(path, [name postfix ext]);
  
  %% Specify experimental variables that are included in the model
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
  cfg.behavConditions = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  % Not sure how to deal with these yet
  cfg.behavConditions = setdiff(cfg.behavConditions, {'trial_nr', 'reward_duration', 'saccade_amplitude'});
  
  %% Support resuming of previously saved fits
  CategorizedModel.nextID(0);
  if lazy && exist(outputFile, 'file')
    newCfg            = cfg;
    load(outputFile);
    if ~isequal(newCfg, keepfield(cfg,fieldnames(newCfg)))
      error('fitSingleCellGLM:cfg', 'Inconsistent analysis configuration found in old model file %s', outputFile);
    end
  else
    categoryModel     = cell(size(data.cellData));
    hierarchicalModel = cell(size(data.cellData));
  end
  
  
  %% Fit model per cell
  parfor iCell = 1:numel(data.cellData)
    if ~isempty(categoryModel{iCell}) && ~isempty(hierarchicalModel{iCell})
      continue
    end
    
    try
      
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                      ...
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)    ...
                      ;
    experiment        = data.experiment;
    experiment.trial  = data.cellData(iCell).trials(selTrials);
    if numel(experiment.trial) < cfg.minNumTrials
      continue;
    end
    
    %% Divide trials into categories according to various stimulus and outcome variables
    trialConditions   = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);
    [trialCategory, ~, categIndex]    = unique(trialConditions, 'rows');
    
    % Construct a joint code for unique data categories
    categDigits       = 1 + max(floor(log10(max(1,trialCategory))), [], 1);
    categCode         = sum(trialCategory .* [10.^cumsum(categDigits(2:end), 'reverse'), 1], 2);
    trialCode         = num2cell(categCode(categIndex));
    [experiment.trial.behavior_code]  = trialCode{:};
    
    %% Specify the form of neural response timecourses to various experimental variables
    modelSpecs        = buildGLM.initDesignSpec(experiment);

    % Responses to behavioral events that vary smoothly with time after the event
    basisFcn          = basisFactory.makeSmoothTemporalBasis('progressive cosine', cfg.eventDuration, cfg.eventNBasis, experiment.binfun);
    for iVar = 1:numel(cfg.behavEvents)
      modelSpecs      = buildGLM.addCovariateTiming(modelSpecs, cfg.behavEvents{iVar}, cfg.behavEvents{iVar}, experiment.desc.(cfg.behavEvents{iVar}), basisFcn);
    end

    %% Specify design matrix using the above parameters
    [design,trialBin] = buildGLM.compileSparseDesignMatrix(modelSpecs, 1:numel(experiment.trial));
    design.dspec.expt.numTimeBins = size(design.X,1);
    design.dspec.expt.trialBins   = trialBin;
    firingRate        = full(buildGLM.getBinnedSpikeTrain(experiment, 'SS', design.trialIndices));
    
    %% Fit model with no specialization of regressors by behavior categories
    fullModel         = CategorizedModel(firingRate, design);
    fullModel.regress(cfg.nCVFolds, cfg.nLambdas);
    
    %% Fit model separately for data per unique category
    categModel        = fullModel.specialize('behavior_code', false, cfg.minNumTrials);
    categModel        = categModel{1};
%     parfor iModel = 1:numel(categModel)
%       categModel(iModel).regress(cfg.nCVFolds, cfg.nLambdas, true);
%     end
    categModel.regress(cfg.nCVFolds, cfg.nLambdas, true);
    categoryModel{iCell}      = categModel;
  
    %% Hierarchical model to find a parsimonious set of trial-category-based specializations for this neuron's response
    hierarchicalModel{iCell}  = buildHierarchicalModel(fullModel, cfg.behavConditions, cfg);
    
    catch err
      displayException(err);
    end
  end
  
  %%
  keyboard
  save(outputFile, 'categoryModel', 'hierarchicalModel', 'cfg');
  
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
    
    %% Parallel optimization for all candidate models
    parfor iCand = 1:numel(candModel)
      candModel{iCand}      = candModel{iCand}.regress(cfg.nCVFolds, cfg.nLambdas, true);
    end
    
    %% Select model with the best relative likelihood w.r.t. the less specialized parent model
    [sortedRelLikeli,iOrder]= sort(cellfun(@(x) x.relativeLikelihood(), candModel), 'ascend');
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
