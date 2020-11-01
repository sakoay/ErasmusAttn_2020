function fitBehaviorDependence(dataFile, postfix)

  %% Default arguments
  if ~exist('postfix', 'var') || isempty(postfix)
    postfix           = '';
  end

  %% Analysis configuration
  cfg                 = struct();
  
  % Cross validation and permutation tests
  cfg.nBootstraps     = 1000;
  cfg.nCVFolds        = 5;
  cfg.nLambdas        = 20;
  
  % Trial selection criteria
  cfg.maxCueStart     = 300;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
  cfg.numHistoryTerms = 2;              % number of past-trial conditions to include in model
  
  % Configuration for lassoglm()
  cfg.target          = 'condition_code';
  cfg.fitOptions      = {'binomial', 'Link', 'logit', 'MaxIter', 1e5, 'Alpha', 0.5};
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('behavModel_%s_%dpast', regexprep(name, '^[^_]+_', ''), cfg.numHistoryTerms);
  outputFile          = fullfile(path, [name postfix ext]);
  
  %% Specify experimental variables that are included in the model
  cfg.stimulusConditions  = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  cfg.stimulusConditions  = setdiff(cfg.stimulusConditions, {'trial_nr', 'condition_code', 'past_condition', 'reward_duration', 'saccade_amplitude'});
  cfg.pastConditions  = [{'condition_code', 'reward_duration'}, cfg.stimulusConditions];

  
  %% Configure models per cell
  experiment          = cell(numel(data.cellData) + 1,1);
  Xvariables          = cell(numel(data.cellData) + 1,1);
  behaviorX           = cell(numel(data.cellData) + 1,1);
  targetY             = cell(numel(data.cellData) + 1,1);
  for iCell = 1:numel(data.cellData)
    %% Select data for a single cell
    experiment{iCell}        = data.experiment;
    experiment{iCell}.id     = [data.cellData(iCell).monkey '_' data.cellData(iCell).cell_id];
    experiment{iCell}.trial  = data.cellData(iCell).trials;
    
    %% Check that trials are temporally ordered
    trialIndex        = [experiment{iCell}.trial.trial_nr];
    if ~issorted(trialIndex, 'strictascend')
      error('fitBehaviorDependence:trial_nr', 'Trials do not have ascending indices.');
    end
    
    %% Define past trial conditions
    trialConditions   = accumfun(2, @(x) cat(1,experiment{iCell}.trial.(x)), cfg.pastConditions);
    pastConditions    = nan(trialIndex(end), numel(cfg.pastConditions), cfg.numHistoryTerms);
    for iPast = 1:cfg.numHistoryTerms
      iTrial          = trialIndex + iPast;
      sel             = iTrial <= trialIndex(end);
      pastConditions(iTrial(sel),:,iPast) = trialConditions(sel,:);
    end
    pastConditions    = pastConditions(trialIndex,:,:);
    
    %% Apply trial selection to cell data, including the requirement that trials have the desired number of historical info
    selTrials         = ismember([experiment{iCell}.trial.condition_code], cfg.selectConditions)    ...
                      & ~any(isnan(pastConditions),2:3)'                                            ...
                      ;
    experiment{iCell}.selectedTrials     = selTrials;
    experiment{iCell}.trial(~selTrials)  = [];
    pastConditions(~selTrials,:,:)= [];
    
    if isempty(experiment{iCell}.trial)
      continue;
    end
    
    targetY{iCell}    = cat(1, experiment{iCell}.trial.(cfg.target));
    
    %% Define regressors; N.B. this assumes that all variables of interest in cfg.stimulusConditions are categorical
    trialConditions   = accumfun(2, @(x) cat(1,experiment{iCell}.trial.(x)), cfg.stimulusConditions);
    [stimCategory,~,stimIndex]    = unique(trialConditions, 'rows');
    stimLabel         = accumfun(2, @(x) strcat(regexprep( cfg.stimulusConditions{x},'_.*',''), '='                                                                 ...
                                                         , parseVariableSpecifications(experiment{iCell}.desc.(cfg.stimulusConditions{x}),stimCategory(:,x),true)   ...
                                                         ), 1:numel(cfg.stimulusConditions));
    
    varLabels         = cell(1, size(stimCategory,1));
    behavX            = nan(size(trialConditions,1), size(stimCategory,1));
    for iStim = 1:size(stimCategory,1)
      varLabels{iStim}= strjoin(stimLabel(iStim,:), ', ');
      behavX(:,iStim) = stimIndex == iStim;
    end
    
    %% Define history-based regressors
    for iCond = 1:numel(cfg.pastConditions)
      varSpecs        = parseVariableSpecifications(experiment{iCell}.desc.(cfg.pastConditions{iCond}));
      [isStim,stimIdx]= ismember(cfg.pastConditions{iCond}, cfg.stimulusConditions);
      if isStim
        %% Convert stimulus specifiers to a test of whether the stimulus has changed
        for iPast = 1:cfg.numHistoryTerms
          varLabels{end+1}        = sprintf('past%d %s same', iPast, regexprep(cfg.pastConditions{iCond},'[_ ].*',''));
          behavX(:,end+1)         = pastConditions(:,iCond,iPast) == trialConditions(:,stimIdx);
        end
      elseif ~isempty(varSpecs)
        %% Treat categorical quantities as indicator variables
        category      = unique(pastConditions(:,iCond,1));
        for iPast = 1:cfg.numHistoryTerms
          for iCat = 1:numel(category)
            varLabels{end+1}      = sprintf('past%d %s', iPast, regexprep(varSpecs{strcmp(varSpecs(:,1),num2str(category(iCat))),2},'[_ ].*',''));
            behavX(:,end+1)       = pastConditions(:,iCond,iPast) == category(iCat);
          end
        end
      else
        %% Keep other regressors as-is
        for iPast = 1:cfg.numHistoryTerms
          varLabels{end+1}        = sprintf('past%d %s', iPast, regexprep(cfg.pastConditions{iCond},'[_ ].*',''));
          behavX(:,end+1)         = pastConditions(:,iCond,iPast);
        end
      end
    end
    
    %% Store regressors 
    behaviorX{iCell}              = behavX;
    Xvariables{iCell}             = varLabels;
%     longfigure; imagesc(zscore(behavX)'); ticklabels(varLabels,'y'); tempscale
  end
  
  %% Pool data across cells
  [nVars,iRef]        = max(cellfun(@numel, Xvariables));
  experiment{end}     = data.experiment.id;
  Xvariables{end}     = Xvariables{iRef};
  behaviorX{end}      = zeros(sum(cellfun(@(x) size(x,1), behaviorX)), nVars);      % N.B. assume that data without these variables have value zero for the corresponding regressor
  range               = 0;
  for iCell = 1:numel(data.cellData)
    if isempty(behaviorX{iCell})
      continue
    end
    range             = range(end) + (1:size(behaviorX{iCell},1));
    behaviorX{end}(range, ismember(Xvariables{end}, Xvariables{iCell}))    ...
                      = behaviorX{iCell};
  end
  targetY{end}        = cat(1, targetY{1:end-1});
  
%   longfigure; imagesc(zscore(behaviorX{end})'); ticklabels(Xvariables{end},'y'); tempscale
  
  %% Logistic regression
  behaviorModel       = cell(numel(data.cellData) + 1,1);
  for iData = 1:numel(behaviorModel)
    %% Require a minimum number of trials for fitting the model
    if numel(targetY{iData}) < cfg.nCVFolds * size(behaviorX{iData},2)
      continue
    end
    
    %% Use cross-validation to determine regularization strength
    X                 = zscore(behaviorX{iData}, 0, 1);
    Y                 = targetY{iData} > 1;
    [modelW,fitInfo]  = lassoglm(X, Y, cfg.fitOptions{:}, 'NumLambda', cfg.nLambdas, 'CV', cfg.nCVFolds);
    lambdaRatio       = fitInfo.LambdaMinDeviance / fitInfo.Lambda(end);
%     lassoPlot(modelW,fitInfo,'plottype','CV');

    %% Use bootstrap experiments to obtain uncertainties on fitted weights
    [~,bootstrapExp]  = permutationExperiments(numel(Y), 0, cfg.nBootstraps, cfg.nCVFolds);
    bootstrapExp(:,1) = [];               % omit original experiment
    if lambdaRatio < 1
      bootstrapW      = nan(size(modelW,1), cfg.nBootstraps);
      parfor iBS = 1:cfg.nBootstraps
        bsW           = lassoglm(X(bootstrapExp(:,iBS),:), Y(bootstrapExp(:,iBS)), cfg.fitOptions{:}, 'NumLambda', 2, 'LambdaRatio', lambdaRatio);
        bootstrapW(:,iBS) = bsW(:,1);
      end
    else
      bootstrapW      = zeros(size(modelW,1), cfg.nBootstraps);
    end
    
    %% Collect fit information
    fitInfo.Xvariables= Xvariables{iData};
    fitInfo.designX   = X;
    fitInfo.targetY   = Y;
    fitInfo.modelW    = modelW;
    fitInfo.bootstrapW= bootstrapW;
    behaviorModel{iData}  = fitInfo;
  end
  
  %% Save output
  save(outputFile, 'experiment', 'behaviorModel', 'cfg');
  fprintf(' -->  %s\n', outputFile);
  
end
