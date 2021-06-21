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
%   cfg.nCVFolds        = 5;
%   cfg.nShuffles       = 100;
  cfg.nCVFolds        = 3;
  cfg.nShuffles       = 50;
  
  % Trial selection criteria
  cfg.maxCueStart     = 1000;            % maximum allowed interval from fixation on to cue presentation, in ms
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
  cfg.selectPastCond  = 1:2;            % keep only trials with these past_condition
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
  
  %% Specify behavioral timings to which the neural data is aligned
  cfg.behavEvents     = rowvec(fieldswithvalue(data.experiment.type, 'timing'));
%   cfg.behavEvents     = setdiff(cfg.behavEvents, {'saccade_onset'});      % remove timings that can have NaNs in particular trials

  cfg.behavEvents(~cellfun(@isempty, regexp(cfg.behavEvents, '_stop$', 'once')))  = [];

  %% Specify experimental variables that are included in the model
  cfg.behavCategories = rowvec(fieldswithvalue(data.experiment.type, 'value'));
%   cfg.behavCategories = setdiff(cfg.behavCategories, {'trial_nr', 'reward_duration', 'saccade_amplitude', 'saccade_direction', 'past_saccade', 'past_condition', 'past_C'});
  cfg.behavCategories = setdiff(cfg.behavCategories, {'trial_nr', 'reward_duration', 'saccade_amplitude', 'condition_code', 'past_condition', 'past_C', 'C_location'});
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('decoding_%s_pastSel_%dvar_%dms_%dshuffle', regexprep(regexprep(name,'^[^_]+_',''),'_[0-9]+ms',''), numel(cfg.behavCategories), cfg.timeBin, cfg.nShuffles);
  if ~isempty(baselineTrials)
    name              = sprintf('%s_baseline%d', name, baselineTrials);
  end
  outputFile          = fullfile(path, [name ext]);
  
  %% ======================================================================
  
  %% Apply trial-level selections
%   longfigure; hold on; plot(arrayfun(@(x) numel(x.trials), data.cellData))
  for iCell = 1:numel(data.cellData)
    data.cellData(iCell).trials( ~ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions) )  = [];
    data.cellData(iCell).trials( ~ismember([data.cellData(iCell).trials.past_condition], cfg.selectPastCond  ) )  = [];
  end
%   plot(arrayfun(@(x) numel(x.trials), data.cellData))
  
  
  %% UGLY : use the entire dataset to determine legal combinations of cfg.behavCategories
  trialCategory       = accumfun(1, @(y) accumfun(2, @(x) cat(1,y.trials.(x)), cfg.behavCategories), data.cellData);
  assert(~any(isnan(trialCategory(:))));
  
  exptCategories      = unique(trialCategory, 'rows');
  
  %% Deduce allowed subsets of task condition values for which we can fit a single decoder
  subConditions       = cell(size(cfg.behavCategories));
  for iCond = 1:numel(cfg.behavCategories)
    %% Get the set of possible values of this behavioral condition, and all possible categories of the other conditions
    [condValues,~,condIdx]    = unique(exptCategories(:,iCond));
    [otherCat,~,otherIndex]   = unique(exptCategories(:,[1:iCond-1,iCond+1:end]), 'rows');
    
    %% Group condition values into sets that all have the same list of other conditions
    allowedCategories         = arrayfun(@(x) otherIndex(condIdx==x), 1:size(condValues,1), 'UniformOutput', false);
    [otherList,iOther,iValue] = unique(cellfun(@(x) mat2str(sort(x)), allowedCategories, 'UniformOutput', false));
    subConditions{iCond}      = arrayfun(@(x) rowvec(condValues(iValue == x)), 1:numel(iOther), 'UniformOutput', false);
    
    %% Must have at least two values to decode
    assert(all( cellfun(@numel, subConditions{iCond}) > 1 ));
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
                      & [data.cellData(iCell).trials.duration] > [data.cellData(iCell).trials.reward_start] + 200   ... HACK for incorrect trial durations
                      & all(accumfun(1, @(x) isfinite([data.cellData(iCell).trials.(x)]), cfg.behavEvents), 1)      ... must have timing info for all events
                      & arrayfun(@(x) ~isempty(x.SS), data.cellData(iCell).trials(:)')                              ... not sure why some trials can have zero spikes
                      ;
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
      trialMeanFR     = arrayfun(@(x) numel(x.SS)*cfg.timeBin/x.duration, experiment.trial(:));
      baselineFR      = movingAverage(trialMeanFR, baselineTrials, 1);
      activityNorm    = 1 ./ baselineFR;
%       longfigure; hold on; yyaxis left; plot(trialMeanFR); yyaxis right; plot(trialMeanFR ./ baselineFR); set(gca,'ylim',[0 2])
%       baselineFR      = movingAverage(trialMeanFR, baselineTrials, 1, [], [], true);
%       plot(trialMeanFR ./ baselineFR,'k-'); set(gca,'ylim',[0 2])
    end
    
    %% Compute binned firing rates aligned to all specified behavioral event times
    alignedActivity   = cell(1, numel(cfg.behavEvents));
    alignedTime       = cell(1, numel(cfg.behavEvents));
    for iEvent = 1:numel(cfg.behavEvents)
      %% Determine bin index for spikes aligned to a particular event time
      rateBins        = arrayfun(@(x) colvec(floor((x.SS - x.(cfg.behavEvents{iEvent})) / cfg.timeBin)), experiment.trial, 'UniformOutput', false);
      minMaxBin       = extrema(cat(1, rateBins{:}));
      binOffset       = 1 - minMaxBin(1);
      
      %% Sum the number of spikes that fall within each time bin, per trial
      alignedTime{iEvent}     = (minMaxBin(1):minMaxBin(2)) * cfg.timeBin;
      nTimeBins               = numel(alignedTime{iEvent});
      alignedActivity{iEvent} = nan(numel(experiment.trial), nTimeBins);
      for iTrial = 1:numel(rateBins)
        alignedActivity{iEvent}(iTrial,:) = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
      end
      
      %% Truncate time bins that don't have data across all trials
      trialStart      = arrayfun(@(x) floor(-x.(cfg.behavEvents{iEvent})                 / cfg.timeBin), experiment.trial);
      trialStop       = arrayfun(@(x) floor((x.duration-1 - x.(cfg.behavEvents{iEvent})) / cfg.timeBin), experiment.trial);
      trialRange      = binOffset + (max(trialStart):min(trialStop));
      
%       longfigure; imagesc(alignedActivity{iEvent}(:,trialRange)'); colorbar
%       longfigure; imagesc(alignedActivity{iEvent}(:,trialRange)' .* activityNorm'); colorbar
%       figure; bandplot(alignedActivity{iEvent}(:,trialRange)'); figure; bandplot(alignedActivity{iEvent}(:,trialRange)' .* activityNorm');
      
      alignedActivity{iEvent} = alignedActivity{iEvent}(:,trialRange) .* activityNorm;
      alignedTime{iEvent}     = alignedTime{iEvent}(trialRange);
    end

    %% Unused because too complicated : Additionally compute firing rates in warped time
%     [alignedActivity{end+1}, alignedTime{end}]  = warpedTimeFiringRate(experiment.trial, cfg.behavEvents, cfg.timeBin);
    
    %% Fit decoder separately per condition type 
    decodeModel       = cell(size(cfg.behavCategories));
    for iCond = 1:numel(cfg.behavCategories)
      %% Restrict decoding to the subset of trials with a design that can be balanced across conditions
      condDecode      = repmat(struct('variable', cfg.behavCategories{iCond}), size(subConditions{iCond}));
      for iSub = 1:numel(subConditions{iCond})
        %% Get list of all combinations of other categories that occur in conjunction with the variable to be decoded
        selTrials     = ismember(trialCategory(:,iCond), subConditions{iCond}{iSub});
        [otherCat,~,otherIndex] = unique(trialCategory(selTrials,[1:iCond-1, iCond+1:end]), 'rows');
        
        %% Count the number of trials in each of the other-category combinations (rows)
        % Columns are for different possible values of the condition to be decoded
        nInOtherCat   = accumfun(2, @(x) sumDataByBin([], otherIndex(trialCategory(selTrials,iCond) == x), size(otherCat,1)), subConditions{iCond}{iSub});
        
        %% Discard trials for which the other-category combination does not occur in conjunction with all values of the variable to be decoded
        discardCat    = any(nInOtherCat < 1, 2);
        doDiscard     = ismember(otherIndex, find(discardCat));
        otherIndex(doDiscard) = [];
        selTrials(selTrials)  = ~doDiscard;
        
        %% Weigh trials so that every category has equal weight, which means that the weighted covariance matrix is the identity matrix
        catWeight     = 1 ./ nInOtherCat;
        trialWeight   = nan(size(otherIndex));
        for iValue = 1:numel(subConditions{iCond}{iSub})
          subSel      = trialCategory(selTrials,iCond) == subConditions{iCond}{iSub}(iValue);
          trialWeight(subSel) = catWeight(otherIndex(subSel), iValue);
        end
        trialWeight   = trialWeight / sum(trialWeight);

        %% Unused : sanity check that there is no remaining correlation of the variable values w.r.t. the other-category combinations
        %{
        [~,~,iOther]  = unique(trialCategory(selTrials,[1:iCond-1, iCond+1:end]), 'rows');
        corr([trialCategory(selTrials,iCond), iOther])
        weightedcorrs([trialCategory(selTrials,iCond), iOther], trialWeight)
        %}
        
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
        parfor iEvent = 1:numel(neuralX)
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
      decodeModel{iCond}              = condDecode;
    end
    
    %% Data format
    decoder{iCell}    = struct('experiment', experiment, 'exptCategories', exptCategories);
    for iCond = 1:numel(cfg.behavCategories)
      decoder{iCell}.(cfg.behavCategories{iCond}) = decodeModel{iCond};
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
