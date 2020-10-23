function fitSingleCellGLM(dataFile)
  
  %% Analysis configuration
  cfg                 = struct();
  
  % Trial selection criteria
  cfg.maxCueStart     = 300;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
  
  % Maximum duration of neural responses per behavioral event
  cfg.eventDuration   = 1500;           % in ms
  cfg.eventNBasis     = 4;              % number of basis functions to use per event-response
  
  % Cross validation and permutation tests
  cfg.nShuffles       = 0;
  cfg.nBootstrap      = 0;
  cfg.nCVFolds        = 10;
  
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
    model             = buildGLM.initDesignSpec(experiment);

    % Responses to behavioral events that vary smoothly with time after the event
    basisFcn          = basisFactory.makeSmoothTemporalBasis('progressive cosine', cfg.eventDuration, cfg.eventNBasis, experiment.binfun);
    for variable = cfg.behavEvents
      model           = buildGLM.addCovariateTiming(model, variable{:}, variable{:}, experiment.desc.(variable{:}), basisFcn);
    end
    
    
    %% Build design matrix 
    design            = buildGLM.compileSparseDesignMatrix(model, 1:numel(experiment.trial), trialConditions);
    design            = buildGLM.removeConstantCols(design);
%     design            = buildGLM.addBiasColumn(design);           % allow systematic offset from zero
    
    firingRate        = full(buildGLM.getBinnedSpikeTrain(experiment, 'SS', design.trialIndices));
    trialLength       = SplitVec(design.trialX, 'equal', 'length');
    
%     longfigure; imagesc(design.X'); longfigure; plot(firingRate)

    %%
    [shuffleExp, bootstrapExp, cvTrainSel] = permutationExperiments(trialLength, cfg.nShuffles, cfg.nBootstrap, cfg.nCVFolds);
    fitInfo           = linearSupportVectorMR( firingRate, design.X, cvTrainSel, shuffleExp, bootstrapExp, [], [], [], true );
  end
  
end
