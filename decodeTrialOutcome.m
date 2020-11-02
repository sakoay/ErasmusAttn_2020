function decodeTrialOutcome(dataFile, postfix)

  %% Default arguments
  if ~exist('postfix', 'var') || isempty(postfix)
    postfix           = '';
  end

  %% Analysis configuration
  cfg                 = struct();
  
  % Cross validation and permutation tests
  cfg.nCVFolds        = 5;
  cfg.nShuffles       = 100;
  
  % Trial selection criteria
  cfg.maxCueStart     = 1000;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
%   cfg.selectPastCond  = 1:2;            % keep only trials with these past_condition
  cfg.minNumTrials    = 20;             % minimum number of trials for fitting models
  
  % Configuration for fitcsvm()
  cfg.timeBin         = 100;            % time binning in ms
  cfg.fitOptions      = {'Verbose', 0, 'Crossval', 'on', 'KFold', cfg.nCVFolds, 'Standardize', true};
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('decoding_%s_%dms_min%dtrials', regexprep(regexprep(name,'^[^_]+_',''),'_[0-9]+ms',''), cfg.timeBin, cfg.minNumTrials);
  outputFile          = fullfile(path, [name postfix ext]);
  
  %% Specify behavioral timings to which the neural data is aligned
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
%   cfg.behavEvents     = setdiff(cfg.behavEvents, {'saccade_onset'});      % remove timings that can have NaNs in particular trials
  
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
  fprintf('Decoding behavioral outcomes from neural activity:\n');
  decoder             = cell(size(data.cellData));
  for iCell = 1:numel(data.cellData)
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                              ...
                      & [data.cellData(iCell).trials.duration] > [data.cellData(iCell).trials.reward_start] + 2*cfg.timeBin   ...
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)            ...
                      & all(accumfun(1, @(x) isfinite([data.cellData(iCell).trials.(x)]), cfg.behavEvents), 1)  ...
                      ;
%                       & ismember([data.cellData(iCell).trials.past_condition], cfg.selectPastCond)      ...
    experiment        = data.experiment;
    experiment.id     = [data.cellData(iCell).monkey '_' data.cellData(iCell).cell_id];
    experiment.trial  = data.cellData(iCell).trials(selTrials);
    experiment.selectedTrials = selTrials;
    
    if numel(experiment.trial) < cfg.minNumTrials
      decoder{iCell}  = experiment;
      continue;
    end
    
    tStart            = tic;
    fprintf('  [%3d]  %-20s (%3d/%-4d = %.3g%% trials)', iCell, experiment.id, sum(selTrials), numel(selTrials), 100*mean(selTrials));
    drawnow;

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
    
    catWeight         = 1 ./ nInCategory;
    trialWeight       = catWeight(catIndex);
    trialWeight       = trialWeight(:) / sum(trialWeight);
    
    %% Compute binned firing rates aligned to all specified behavioral event times
    alignedActivity   = cell(size(cfg.behavEvents));
    alignedTime       = cell(size(cfg.behavEvents));
    for iEvent = 1:numel(cfg.behavEvents)
      trials          = experiment.trial;
      
      %% Determine bin index for spikes aligned to a particular event time
      rateBins        = arrayfun(@(x) floor((x.SS - x.(cfg.behavEvents{iEvent})) / cfg.timeBin), trials, 'UniformOutput', false);
      minMaxBin       = extrema(cat(1,rateBins{:}));
      binOffset       = 1 - minMaxBin(1);
      
      %% Sum the number of spikes that fall within each time bin, per trial
      alignedTime{iEvent}     = (minMaxBin(1):minMaxBin(2)) * cfg.timeBin;
      nTimeBins               = numel(alignedTime{iEvent});
      alignedActivity{iEvent} = nan(numel(rateBins), nTimeBins);
      for iTrial = 1:numel(rateBins)
        alignedActivity{iEvent}(iTrial,:) = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
      end
      
      %% Truncate time bins that don't have data across all trials
      trialStart      = arrayfun(@(x) floor(-x.(cfg.behavEvents{iEvent})                 / cfg.timeBin), trials);
      trialStop       = arrayfun(@(x) floor((x.duration-1 - x.(cfg.behavEvents{iEvent})) / cfg.timeBin), trials);
      trialRange      = binOffset + (max(trialStart):min(trialStop));
      alignedActivity{iEvent} = alignedActivity{iEvent}(:,trialRange);
      alignedTime{iEvent}     = alignedTime{iEvent}(trialRange);
%       figure; imagesc(alignedActivity{iEvent})
    end

    %% Define permutation experiments
    shuffleExp        = permutationExperiments(numel(experiment.trial), cfg.nShuffles, [], []);
    decoder{iCell}    = struct('experiment', experiment, 'shuffleExp', shuffleExp);
    
    %% Fit decoder separately per condition type 
    for iCond = 1:numel(cfg.behavConditions)
      behaviorY       = trialCondition(:,iCond);
      
      %% Fit decoder for neural data aligned to different behavioral events
      svModel         = cell(1, numel(cfg.behavEvents));
      parfor iEvent = 1:numel(cfg.behavEvents)
        %% Fit decoder per time bin
        nTimeBins     = size(alignedActivity{iEvent}, 2);
        model         = struct('prediction', nan([size(alignedActivity{iEvent}), size(shuffleExp,2)]));
        for iExp = 1:size(shuffleExp,2)
          sel         = shuffleExp(:,iExp);
          for iTime = 1:nTimeBins
            fitInfo   = fitcsvm(alignedActivity{iEvent}(sel,iTime), behaviorY, cfg.fitOptions{:}, 'Weights', trialWeight);
            model.prediction(:,iTime,iExp)  = fitInfo.kfoldPredict();
          end
        end
        
        %% Define classification accuracy using weighted trials (kfoldLoss() does not seem to account for this)
        model.accuracy= sqsum(trialWeight .* (model.prediction == behaviorY), 1);
        model.prediction(:,:,2:end) = [];       % reduce data storage
        svModel{iEvent} = model;
      end
      svModel         = [svModel{:}];
      
      %% Collect fit information
      [svModel.alignment]       = cfg.behavEvents{:};
      [svModel.alignedTime]     = alignedTime{:};
      [svModel.alignedActivity] = alignedActivity{:};
      condDecode      = struct('target', behaviorY, 'model', svModel);
      decoder{iCell}.(cfg.behavConditions{iCond}) = condDecode;
    end
    
    %% Progressive save 
%     for ii=1:numel(condDecode.model); mm=condDecode.model(ii); figure; bandplot(mm.alignedTime,mm.accuracy(:,2:end)); hold on; plot(mm.alignedTime,mm.accuracy(:,1),'k'); xlabel(mm.alignment); end; tilefigs
    
    save(outputFile, 'decoder', 'cfg');
    fprintf(' ... %gs\n', toc(tStart));
  end
  
  
  %% Save models
  fprintf(' -->  %s\n', outputFile);
  
end
