function decodeTrialOutcome(dataFile, postfix)

  %% Default arguments
  if ~exist('postfix', 'var') || isempty(postfix)
    postfix           = '';
  end

  %% Analysis configuration
  cfg                 = struct();
  
  % Cross validation and permutation tests
  cfg.nCVFolds        = 5;
  
  % Trial selection criteria
  cfg.maxCueStart     = 300;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
%   cfg.selectPastCond  = 1:2;            % keep only trials with these past_condition
  cfg.minNumTrials    = 20;             % minimum number of trials for fitting models
  
  % Configuration for fitcsvm()
  cfg.fitOptions      = {'Verbose', 0, 'Crossval', 'on', 'KFold', cfg.nCVFolds};
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('outcomeDecode_%s_min%dtrials', regexprep(name, '^[^_]+_', ''), cfg.minNumTrials);
  outputFile          = fullfile(path, [name postfix ext]);
  
  %% Specify behavioral timings to which the neural data is aligned
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
  
  %% Specify experimental variables that are included in the model
  cfg.behavCategories = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  cfg.behavCategories = setdiff(cfg.behavCategories, {'trial_nr', 'reward_duration', 'saccade_amplitude', 'past_condition'});
  sel                 = contains(cfg.behavCategories, 'condition');
  cfg.behavConditions = cfg.behavCategories(sel);
  cfg.behavCategories(sel)  = [];
  
  %% UGLY : use the entire dataset to determine legal combinations of cfg.behavCategories
  trialCategory       = accumfun(1, @(y) accumfun(2, @(x) cat(1,y.trials.(x)), [cfg.behavCategories,cfg.behavConditions]), data.cellData);
  exptCategories      = unique(trialCategory, 'rows');
  
  % Apply the specified behavior-level selections 
  [~,iCond]           = ismember(cfg.behavConditions, 'condition_code');
  if ~isempty(iCond)
    exptCategories(~ismember(exptCategories(:,numel(cfg.behavCategories) + iCond), cfg.selectConditions), :)  = [];
  end
  
  
  %% Fit decoding models per cell
  decoder             = cell(size(data.cellData));
  for iCell = 1:numel(data.cellData)
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                      ...
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)    ...
                      ;
%                       & ismember([data.cellData(iCell).trials.past_condition], cfg.selectPastCond)      ...
    experiment        = data.experiment;
    experiment.id     = [data.cellData(iCell).monkey '_' data.cellData(iCell).cell_id];
    experiment.trial  = data.cellData(iCell).trials(selTrials);
%     mean(selTrials)
    
    if numel(experiment.trial) < cfg.minNumTrials
      decoder{iCell}  = experiment;
      continue;
    end
    
    %% Define behavioral variables to decode from neural activity
    trialCategory     = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavCategories);
    trialCondition    = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);

    %% Weigh trials so that every category has equal weight, which means that the weighted covariance matrix is the identity matrix
    [~,catIndex]      = ismember([trialCategory, trialCondition], exptCategories, 'rows');
    nInCategory       = arrayfun(@(x) sum(catIndex==x), 1:size(exptCategories,1));
    if any(nInCategory < 1)
      warning('decodeTrialOutcome:reweight', 'Insufficient trials to have data in each experimental category for cell %s, skipping this cell.', experiment.id);
      decoder{iCell}  = experiment;
      continue;
    end
    
    catWeight         = nInCategory / sum(nInCategory);
    trialWeight       = catWeight(catIndex);
    
    %% Compute binned firing rates aligned to all specified behavioral event times
    selTrials         = cell(size(cfg.behavEvents));
    alignedActivity   = cell(size(cfg.behavEvents));
    alignedTime       = cell(size(cfg.behavEvents));
    for iEvent = 1:numel(cfg.behavEvents)
      %% Include only those trials with valid event times of this type
      selTrials{iEvent} = arrayfun(@(x) isfinite(x.(cfg.behavEvents{iEvent})), experiment.trial);
      trials          = experiment.trial(selTrials{iEvent});
      
      %% Determine bin index for spikes aligned to a particular event time
      rateBins        = arrayfun(@(x) floor((x.SS - x.(cfg.behavEvents{iEvent})) / experiment.binSize), trials, 'UniformOutput', false);
      minMaxBin       = extrema(cat(1,rateBins{:}));
      binOffset       = 1 - minMaxBin(1);
      
      %% Sum the number of spikes that fall within each time bin, per trial
      alignedTime{iEvent}     = (minMaxBin(1):minMaxBin(2)) * experiment.binSize;
      nTimeBins               = numel(alignedTime{iEvent});
      alignedActivity{iEvent} = nan(numel(rateBins), nTimeBins);
      for iTrial = 1:numel(rateBins)
        alignedActivity{iEvent}(iTrial,:) = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
      end
      
      %% Truncate time bins that don't have data across all trials
      trialStart      = arrayfun(@(x) floor(-x.(cfg.behavEvents{iEvent})                 / experiment.binSize), trials);
      trialStop       = arrayfun(@(x) floor((x.duration-1 - x.(cfg.behavEvents{iEvent})) / experiment.binSize), trials);
      trialRange      = binOffset + (max(trialStart):min(trialStop));
      alignedActivity{iEvent} = alignedActivity{iEvent}(:,trialRange);
      alignedTime{iEvent}     = alignedTime{iEvent}(trialRange);
%       figure; imagesc(alignedActivity{iEvent})
    end
    
    %% Fit decoder separately per condition type 
    for iCond = 1:numel(cfg.behavConditions)
      behaviorY       = trialCondition(:,iCond);
      
      %% Fit decoder for neural data aligned to different behavioral events
      tic
      svModel         = cell(1, numel(cfg.behavEvents));
      parfor iEvent = 1:numel(cfg.behavEvents)
        %% Apply event-specific trial selection
        behavY        = behaviorY(selTrials{iEvent});
        trialW        = trialWeight(selTrials{iEvent});
        
        %% Fit decoder per time bin
        nTimeBins     = size(alignedActivity{iEvent}, 2);
        fit           = cell(nTimeBins, 1);
        for iTime = 1:nTimeBins
          fit{iTime}  = fitcsvm(alignedActivity{iEvent}(:,iTime), behavY, cfg.fitOptions{:}, 'Weights', trialW);
        end
        svModel{iEvent} = fit;
      end
      toc
    end
  end
  
  
  %% Save models
  save(outputFile, 'unspecializedModel', 'categoryModel', 'hierarchicalModel', 'cfg');
  fprintf(' -->  %s\n', outputFile);
  
end
