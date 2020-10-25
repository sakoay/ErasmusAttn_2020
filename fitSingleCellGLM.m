function fitSingleCellGLM(dataFile)
  
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
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Specify experimental variables that are included in the model
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
  cfg.behavConditions = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  % Not sure how to deal with these yet
  cfg.behavConditions = setdiff(cfg.behavConditions, {'trial_nr', 'reward_duration', 'saccade_amplitude'});
    
  %% Fit model per cell
  for iCell = 1:numel(data.cellData)
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                      ...
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)    ...
                      ;
    experiment        = data.experiment;
    experiment.trial  = data.cellData(iCell).trials(selTrials);
    
    %% Divide trials into categories according to various stimulus and outcome variables
    trialConditions   = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);
    [trialCategory, ~, categIndex]  = unique(trialConditions, 'rows');
    
    %% Specify the form of neural response timecourses to various experimental variables
    modelSpecs        = buildGLM.initDesignSpec(experiment);

    % Responses to behavioral events that vary smoothly with time after the event
    basisFcn          = basisFactory.makeSmoothTemporalBasis('progressive cosine', cfg.eventDuration, cfg.eventNBasis, experiment.binfun);
    for variable = cfg.behavEvents
      modelSpecs      = buildGLM.addCovariateTiming(modelSpecs, variable{:}, variable{:}, experiment.desc.(variable{:}), basisFcn);
    end
    
    %% Fit model per trial category
    %{
    fitInfo           = cell(1, size(trialCategory,1));
    parfor iCat = 1:size(trialCategory,1)
      %% Build design matrix 
      design          = buildGLM.compileSparseDesignMatrix(modelSpecs, find(categIndex == iCat));
      design          = buildGLM.removeConstantCols(design);
%     design          = buildGLM.addBiasColumn(design);           % allow systematic offset from zero

      firingRate      = full(buildGLM.getBinnedSpikeTrain(experiment, 'SS', design.trialIndices));
      trialLength     = SplitVec(design.trialX, 'equal', 'length');

%     longfigure; imagesc(design.X'); longfigure; plot(firingRate)

      %% Cross-validated model fit
      [shuffleExp, bootstrapExp, cvTrainSel]  = permutationExperiments(trialLength, cfg.nShuffles, cfg.nBootstrap, cfg.nCVFolds);
      fitInfo{iCat}   = linearSupportVectorMR( firingRate, design.X, cvTrainSel, shuffleExp, bootstrapExp, [], [], [], true );
      fitInfo{iCat}.model = buildGLM.combineWeights(design, fitInfo{iCat}.modelW);
    end
    fitInfo           = [fitInfo{:}];
    
%     sfigure; scatter([fitInfo.shuffleSigni], arrayfun(@(x) x.correlation(1), fitInfo)); xlabel('Significance'); ylabel('corr(model, data)');
    %}
  
    %% Hierarchical model to find a parsimonious set of trial-category-based specializations for this neuron's response
    design            = buildGLM.compileSparseDesignMatrix(modelSpecs, 1:numel(experiment.trial));
    firingRate        = full(buildGLM.getBinnedSpikeTrain(experiment, 'SS', design.trialIndices));
    fullModel         = CategorizedModel(firingRate, design);
    fullModel         = fullModel.regress(cfg.nShuffles, cfg.nCVFolds);
    categorizedModel  = buildHierarchicalModel(fullModel, cfg.behavConditions, cfg);
  end
  
end

function model = buildHierarchicalModel(model, conditionLabels, cfg)
  
  %% Continue specializing models until there is no improvement in goodness of fit
  while ~isempty(conditionLabels)
    %% Generate further specialized models per candidate condition 
    candModel               = cell(size(conditionLabels));
    for iCond = 1:numel(conditionLabels)
      %% Construct models for all possible partitions of values for a given condition label
      tryModels             = model.specialize(conditionLabels{iCond});
      for iTry = numel(tryModels):-1:1
        if any(tryModels{iTry}.numberOfTrials() < cfg.minNumTrials)
          tryModels(iTry)   = [];
        else
          tryModels{iTry}   = tryModels{iTry}.regress(0, cfg.nCVFolds);
        end
      end
      
      %% Select model with the largest effect size 
      [~,iBest]             = max(cellfun(@(x) x.localEffectSize(), tryModels));
      candModel{iCond}      = tryModels{iBest};
    end
    
    %% Check if we should continue specializing models
    candModel(cellfun(@isempty, candModel)) = [];
    if isempty(candModel)
      break;
    end
    [maxEffect,iBest]       = max(cellfun(@(x) x.localEffectSize(), candModel));
    
    %%
    [pan,shape]             = makePanels(numel(fieldnames(model.regressors)));
    iPlot                   = 0;
    for what = fieldnames(model.regressors)'
      iPlot                 = iPlot + 1;
      axs                   = selectPanel(pan, iPlot, shape);
      hold(axs, 'on');
      xlabel(axs, sprintf('time in trial (%s)', model.design.dspec.expt.unitOfTime));
      ylabel(axs, strrep(what{:}, '_', ' '));
      
      for iModel = 1:numel(candModel{iBest})
        plot(axs, candModel{iBest}(iModel).regressors.(what{:}).time, candModel{iBest}(iModel).regressors.(what{:}).response);
      end
      plot(axs, model.regressors.(what{:}).time, model.regressors.(what{:}).response, 'k');
    end
  end
  
end
