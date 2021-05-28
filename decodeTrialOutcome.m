% ff = rdir('C:\Neuroscience\ErasmusAttn\glmExpt_*.mat');
% for ii = 1:numel(ff); decodeTrialOutcome(ff(ii).name); end
function decodeTrialOutcome(dataFile, baselineTrials, lazy)

  %% Default arguments
  if ~exist('baselineTrials', 'var')
    baselineTrials    = 5;
  end
  if ~exist('lazy', 'var') || isempty(lazy)
    lazy              = true;
  end

  %% Analysis configuration
  cfg                 = struct();
  
  % Cross validation and permutation tests
  cfg.nCVFolds        = 5;
  cfg.nShuffles       = 100;
%   cfg.nShuffles       = 50;
  
  % Trial selection criteria
  cfg.maxCueStart     = 1000;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
%   cfg.selectPastCond  = 1:2;            % keep only trials with these past_condition
  cfg.minNumTrials    = 10;             % minimum number of trials for fitting models
  
  % Configuration for model fitting
%   cfg.timeBin         = 150;            % time binning in ms
  cfg.timeBin         = 200;            % time binning in ms
  cfg.codingDesign    = 'onevsone';
  cfg.fitFcn          = { @(X,y,varargin) fitcsvm(X,y, 'Verbose', 0, 'Crossval', 'on', 'KFold', cfg.nCVFolds, 'Standardize', true, varargin{:})                                                       ...
                        , @(X,y,varargin) fitcecoc(X,y, 'Coding', cfg.codingDesign, 'Learners', templateSVM('Standardize',true), 'Verbose', 0, 'Crossval', 'on', 'KFold', cfg.nCVFolds, varargin{:})  ...
                        };
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('decoding_%s_%dms_%dshuffle', regexprep(regexprep(name,'^[^_]+_',''),'_[0-9]+ms',''), cfg.timeBin, cfg.nShuffles);
  if ~isempty(baselineTrials)
    name              = sprintf('%s_baseline%d', name, baselineTrials);
  end
  outputFile          = fullfile(path, [name ext]);
  
  %% Specify behavioral timings to which the neural data is aligned
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
%   cfg.behavEvents     = setdiff(cfg.behavEvents, {'saccade_onset'});      % remove timings that can have NaNs in particular trials

  cfg.behavEvents(~cellfun(@isempty, regexp(cfg.behavEvents, '_stop$', 'once')))  = [];

  %% Specify experimental variables that are included in the model
  cfg.behavCategories = rowvec(fieldswithvalue(data.experiment.type, 'value'));
  cfg.behavCategories = setdiff(cfg.behavCategories, {'trial_nr', 'reward_duration', 'saccade_amplitude', 'past_condition'});
  
  %% UGLY : use the entire dataset to determine legal combinations of cfg.behavCategories
  trialCategory       = accumfun(1, @(y) accumfun(2, @(x) cat(1,y.trials.(x)), cfg.behavCategories), data.cellData);
  exptCategories      = unique(trialCategory, 'rows');
  
  % Apply the specified behavior-level selections 
  [~,iCond]           = ismember('condition_code', cfg.behavCategories);
  exptCategories(~ismember(exptCategories(:,iCond), cfg.selectConditions), :) = [];
  
  %% Deduce allowed subsets of task condition values for which we can fit a single decoder
  subConditions       = cell(size(cfg.behavCategories));
  for iCond = 1:numel(cfg.behavCategories)
    %% Get the set of possible values of this behavioral condition, and all possible categories of the other conditions
    [condValues,~,condIdx]    = unique(exptCategories(:,iCond));
    [otherCat,~,otherIndex]   = unique(exptCategories(:,[1:iCond-1,iCond+1:end]), 'rows');
    
    %% Group condition values by common restrictions on the other categories
    allowedCategories         = arrayfun(@(x) otherIndex(condIdx==x), 1:size(condValues,1), 'UniformOutput', false);
    [~,iOther,iValue]         = unique(cellfun(@(x) mat2str(sort(x)), allowedCategories, 'UniformOutput', false));
    subConditions{iCond}      = arrayfun(@(x) rowvec(condValues(iValue == x)), 1:numel(iOther), 'UniformOutput', false);
  end
  
  
  %% Support resuming of previously saved fits
  if lazy && exist(outputFile, 'file')
    newCfg            = cfg;
    load(outputFile);
    if ~isequal(rmfield(newCfg,'fitFcn'), rmfield(keepfield(cfg,fieldnames(newCfg)),'fitFcn'))
      error('decodeTrialOutcome:cfg', 'Inconsistent analysis configuration found in old model file %s', outputFile);
    end
  else
    decoder           = cell(size(data.cellData));
  end

  %% Fit decoding models per cell
  fprintf('Decoding behavioral outcomes from neural activity:\n');
  for iCell = 1:numel(data.cellData)
    if ~isempty(decoder{iCell})
      continue
    end
    
    %% Apply trial selection to cell data
    selTrials         = [data.cellData(iCell).trials.cue_start] <= cfg.maxCueStart                                  ...
                      & [data.cellData(iCell).trials.duration] > [data.cellData(iCell).trials.reward_start] + 200   ... % HACK for incorrect trial durations
                      & ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions)                ...
                      & all(accumfun(1, @(x) isfinite([data.cellData(iCell).trials.(x)]), cfg.behavEvents), 1)      ...
                      & arrayfun(@(x) ~isempty(x.SS), data.cellData(iCell).trials(:)')                              ... % not sure why some trials can have zero spikes
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

    %% Optionally normalize the firing rates of each trial to the moving average rate 
    if isempty(baselineTrials)
      activityNorm    = 1;
    else
      %%
      trialMeanFR     = arrayfun(@(x) numel(x.SS)/x.duration, experiment.trial(:));
      baselineFR      = movingAverage(trialMeanFR, baselineTrials, 1);
      activityNorm    = 1 ./ baselineFR;
%       longfigure; hold on; yyaxis left; plot(trialMeanFR); yyaxis right; plot(trialMeanFR ./ baselineFR); set(gca,'ylim',[0 2])
%       baselineFR      = movingAverage(trialMeanFR, baselineTrials, 1, [], [], true);
%       plot(trialMeanFR ./ baselineFR,'k-'); set(gca,'ylim',[0 2])
    end
    
    %% Compute binned firing rates aligned to all specified behavioral event times
    alignedActivity   = cell(size(cfg.behavEvents));
    alignedTime       = cell(size(cfg.behavEvents));
    for iEvent = 1:numel(cfg.behavEvents)
      trials          = experiment.trial;
      
      %% Determine bin index for spikes aligned to a particular event time
      rateBins        = arrayfun(@(x) colvec(floor((x.SS - x.(cfg.behavEvents{iEvent})) / cfg.timeBin)), trials, 'UniformOutput', false);
      minMaxBin       = extrema(cat(1, rateBins{:}));
      binOffset       = 1 - minMaxBin(1);
      
      %% Sum the number of spikes that fall within each time bin, per trial
      alignedTime{iEvent}     = (minMaxBin(1):minMaxBin(2)) * cfg.timeBin;
      nTimeBins               = numel(alignedTime{iEvent});
      alignedActivity{iEvent} = nan(numel(trials), nTimeBins);
      for iTrial = 1:numel(rateBins)
        alignedActivity{iEvent}(iTrial,:) = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
      end
      
      %% Truncate time bins that don't have data across all trials
      trialStart      = arrayfun(@(x) floor(-x.(cfg.behavEvents{iEvent})                 / cfg.timeBin), trials);
      trialStop       = arrayfun(@(x) floor((x.duration-1 - x.(cfg.behavEvents{iEvent})) / cfg.timeBin), trials);
      trialRange      = binOffset + (max(trialStart):min(trialStop));
      
%       longfigure; imagesc(alignedActivity{iEvent}(:,trialRange)'); colorbar
%       longfigure; imagesc(alignedActivity{iEvent}(:,trialRange)' .* activityNorm'); colorbar
%       figure; bandplot(alignedActivity{iEvent}(:,trialRange)'); figure; bandplot(alignedActivity{iEvent}(:,trialRange)' .* activityNorm');
      
      alignedActivity{iEvent} = alignedActivity{iEvent}(:,trialRange) .* activityNorm;
      alignedTime{iEvent}     = alignedTime{iEvent}(trialRange);
    end

    
    %% Fit decoder separately per condition type 
    decoder{iCell}    = struct('experiment', experiment, 'exptCategories', exptCategories);
    for iCond = 1:numel(cfg.behavCategories)
      %% Restrict decoding to the subset of trials with a design that can be balanced across conditions
      condDecode      = repmat(struct('variable', cfg.behavCategories{iCond}), size(subConditions{iCond}));
      parfor iSub = 1:numel(subConditions{iCond})
        %% Discard trials if necessary in order to be able to balance trial conditions
        selTrials     = ismember(trialCategory(:,iCond), subConditions{iCond}{iSub});
        condIndex     = [];
        nInCategory   = [];
        for iDiscard = 0:1
          %% Determine all possible combinations of all task condition values for the target subset of trials 
          [otherCat,~,otherIndex] = unique(trialCategory(selTrials,[1:iCond-1, iCond+1:end]), 'rows');
          condValues  = allCombinations(subConditions{iCond}{iSub}, 1:size(otherCat,1));

          %% If any category has no entries, discard trials with the problematic condition values
          [~,condIndex] = ismember([trialCategory(selTrials,iCond), otherIndex], condValues, 'rows');
          nInCategory   = arrayfun(@(x) sum(condIndex==x), 1:size(condValues,1));
          missingValues = condValues(nInCategory < 1, 2);
%           otherCat(condValues(missingValues,2),:)

          %% Apply the additional trial selection
          if isempty(missingValues)
            break
          elseif iDiscard 
            error('decodeTrialOutcome:iDiscard', 'Not sure what to do.');
          else
            selTrials(selTrials)  = ~ismember(otherIndex, missingValues);
          end
        end
        
        
        %% Weigh trials so that every category has equal weight, which means that the weighted covariance matrix is the identity matrix
        catWeight     = 1 ./ nInCategory;
        trialWeight   = catWeight(condIndex);
        trialWeight   = trialWeight(:) / sum(trialWeight);

        %% Define neural and behavioral data
        behaviorY     = trialCategory(selTrials,iCond);
        neuralX       = cellfun(@(x) x(selTrials,:), alignedActivity, 'UniformOutput', false);
        
        condDecode(iSub).values           = unique(behaviorY);
        condDecode(iSub).selectedTrials   = selTrials;
        condDecode(iSub).trialWeight      = trialWeight;
        condDecode(iSub).target           = behaviorY;
        
        if sum(selTrials) < cfg.minNumTrials
          condDecode(iSub).model          = [];
          warning( 'decodeTrialOutcome:selTrials', 'Insufficient trials to have data in each category of %s = %s for cell %s.'  ...
                 , cfg.behavCategories{iCond}, mat2str(condDecode(iSub).values), experiment.id );
          continue;
        end
        
        %% Define permutation experiments
        shuffleExp    = permutationExperiments(numel(behaviorY), cfg.nShuffles, [], []);
        
        %% Fit decoder for neural data aligned to different behavioral events
%         try
%           
%         warning('error', 'stats:cvpartition:KFoldMissingGrp');
        
        optimizer     = cfg.fitFcn{1 + (numel(condDecode(iSub).values) > 2)};
        svModel       = cell(1, numel(cfg.behavEvents));
        for iEvent = 1:numel(cfg.behavEvents)
          %% Fit decoder per time bin
          nTimeBins   = size(neuralX{iEvent}, 2);
          prediction  = nan([size(neuralX{iEvent}), size(shuffleExp,2)]);
          for iExp = 1:size(shuffleExp,2)
            X         = neuralX{iEvent}(shuffleExp(:,iExp),:);
            for iTime = 1:nTimeBins
              fitInfo = optimizer(X(:,iTime), behaviorY, 'Weights', trialWeight);
              prediction(:,iTime,iExp)  = fitInfo.kfoldPredict();
            end
          end

          %% Define classification accuracy using weighted trials (kfoldLoss() does not seem to account for this)
          model       = struct('prediction', prediction);
          model.accuracy              = sqsum(trialWeight .* (model.prediction == behaviorY), 1);
          model.prediction(:,:,2:end) = [];       % reduce data storage
          svModel{iEvent}             = model;

%         pp=mean(model.accuracy >= model.accuracy(:,1), 2); figure; plot(pp); hold on; plot(find(pp<=0.05),pp(pp<=0.05),'o'); title(cfg.codingDesign);
%         figure; bandplot(model.accuracy(:,2:end)); hold on; plot(model.accuracy(:,1)); title(cfg.codingDesign);
%         figure; confusionchart(behaviorY,prediction(:,6,1)); title(cfg.codingDesign);
        end
        svModel         = [svModel{:}];
%         warning('on', 'stats:cvpartition:KFoldMissingGrp');

%         catch err
%           keyboard
%         end

        %% Collect fit information
        [svModel.alignment]           = cfg.behavEvents{:};
        [svModel.alignedTime]         = alignedTime{:};
        [svModel.alignedActivity]     = neuralX{:};
        condDecode(iSub).model        = svModel;
      end
      decoder{iCell}.(cfg.behavCategories{iCond}) = condDecode;
    end
    
    %% Progressive save 
%     for ii=1:numel(condDecode.model); mm=condDecode.model(ii); figure; bandplot(mm.alignedTime,mm.accuracy(:,2:end)); hold on; plot(mm.alignedTime,mm.accuracy(:,1),'k'); xlabel(strrep(mm.alignment,'_',' ')); end; tilefigs
    
    save(outputFile, 'decoder', 'cfg');
    fprintf(' ... %gs\n', toc(tStart));
  end
  
  %% Save models
  save(outputFile, 'decoder', 'cfg');
  fprintf(' -->  %s\n', outputFile);
  
end
