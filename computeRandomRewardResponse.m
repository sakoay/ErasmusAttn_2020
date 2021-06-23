% ff = rdir('C:\Neuroscience\ErasmusAttn\glmExpt_*.mat');
% for ii = 1:numel(ff); computeRandomRewardResponse(ff(ii).name); end
function computeRandomRewardResponse(dataFile)

  %% Analysis configuration
  cfg                 = struct();
  
  % Permutation tests
  cfg.nShuffles       = 1000;
  
  % Trial selection criteria
  cfg.minNumTrials    = 10;
  cfg.selectConditions= 4;              % keep only trials with these condition_code
  
  % Configuration for model fitting
  cfg.timeBin         = 100;             % time binning in ms
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('randReward_%s_%dms_%dshuffle', regexprep(regexprep(name,'^[^_]+_',''),'_[0-9]+ms',''), cfg.timeBin, cfg.nShuffles);
  outputFile          = fullfile(path, [name ext]);
  
  %% ======================================================================
  
  
  %% Random reward response quantification per cell
  for iCell = 1:numel(data.cellData)
    %% Apply trial selection to cell data
    selTrials         = ismember([data.cellData(iCell).trials.condition_code], cfg.selectConditions);
    experiment        = data.experiment;
    experiment.id     = [data.cellData(iCell).monkey '_' data.cellData(iCell).cell_id];
    experiment.trial  = data.cellData(iCell).trials(selTrials);
    experiment.selectedTrials = selTrials;
    
    if numel(experiment.trial) < cfg.minNumTrials
      continue;
    end
    
    fprintf('  [%3d]  %-20s (%3d/%-4d = %.3g%% trials)', iCell, experiment.id, sum(selTrials), numel(selTrials), 100*mean(selTrials));
    drawnow;

    %% Compute binned firing rate for experiment 
    spikeTimes        = {experiment.trial.SS};

    % Determine bin index for spikes aligned to the reward time
    rateBins          = cellfun(@(x,y) colvec(floor((x - y) / cfg.timeBin)), spikeTimes, {experiment.trial.reward_start}, 'UniformOutput', false);
    minMaxBin         = extrema(cat(1, rateBins{:}));
    binOffset         = 1 - minMaxBin(1);
    alignedTime       = (minMaxBin(1):minMaxBin(2)) * cfg.timeBin;
    nTimeBins         = numel(alignedTime);
    alignedActivity   = nan(numel(experiment.trial), nTimeBins);    %, 1 + cfg.nShuffles);

    % Sum the number of spikes that fall within each time bin, per trial
    for iTrial = 1:numel(rateBins)
      alignedActivity(iTrial,:,1) = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
    end
    
    %% Truncate time bins that don't have data across all trials
    trialStart        = arrayfun(@(x) floor(-x.reward_start                 / cfg.timeBin), experiment.trial) + 1;
    trialStop         = arrayfun(@(x) floor((x.duration-1 - x.reward_start) / cfg.timeBin), experiment.trial) - 1;
    trialRange        = binOffset + (max(trialStart):min(trialStop));

%       longfigure; imagesc(alignedActivity(:,trialRange)); colorbar

    alignedActivity   = alignedActivity(:,trialRange,:);
    alignedTime       = alignedTime(trialRange);
    
    %%
    trainActivity     = alignedActivity(1:2:end-1,:);
    testActivity      = alignedActivity(2:2:end,  :);
%     stdActivity       = std(trainActivity - mean(trainActivity,2), 0, 'all');
    baseline          = mean( trainActivity(:,alignedTime <  0)           , 2 );
    trainMean         = mean( trainActivity(:,alignedTime >= 0) - baseline, 2 );
    stdActivity       =  std( trainActivity(:,alignedTime <  0) - baseline, 0, 'all' );
    
    rewardOffset      = sum(alignedTime < 0);
    [~,maxTime]       = max(trainMean);         maxTime = rewardOffset + maxTime;
    [~,minTime]       = min(trainMean);         minTime = rewardOffset + minTime;
    
    baseline          = mean( testActivity(:,alignedTime <  0), 2 );
%     stdActivity       =  std( testActivity(:,alignedTime <  0) - baseline, 0, 'all' );
    
    figure; plot(alignedTime, alignedActivity'); xlabel('Time from random reward (ms)'); ylabel('Number of spikes / 100ms');
    figure; plot(alignedTime, (testActivity - baseline)' / stdActivity)
    figure; bandplot(alignedTime,alignedActivity(:,:,1)');  xlabel('Time from random reward (ms)'); ylabel('Number of spikes / 100ms');
    
    %%
    [(testActivity(:,maxTime) - baseline) >  2 * stdActivity, (testActivity(:,minTime) - baseline) < -2 * stdActivity]
    
    %% 
    continue
    
    %% The same but for permutation tests where the spike times are randomly rotated in each trial
    parfor iShuffle = 1:cfg.nShuffles
      %% Randomly rotate spike times for permutation tests
      spikeTimes      = {experiment.trial.SS};
      for iTrial = 1:numel(spikeTimes)
        spikeTimes{iTrial}              = spikeTimes{iTrial} - rand()*experiment.trial(iTrial).duration;
        wraparound                      = spikeTimes{iTrial} < 0;
        spikeTimes{iTrial}(wraparound)  = spikeTimes{iTrial}(wraparound) + experiment.trial(iTrial).duration;
      end
   
      %% Determine bin index for spikes aligned to the reward time
      rateBins        = cellfun(@(x,y) colvec(floor((x - y) / cfg.timeBin)), spikeTimes, {experiment.trial.reward_start}, 'UniformOutput', false);

      %% Sum the number of spikes that fall within each time bin, per trial
      activity        = nan(numel(rateBins),nTimeBins);
      for iTrial = 1:numel(rateBins)
        activity(iTrial,:)              = sumDataByBin([], rateBins{iTrial} + binOffset, nTimeBins);
      end
      alignedActivity(:,:,iShuffle+1)   = activity;
    end
    
    %% Truncate time bins that don't have data across all trials
    trialStart        = arrayfun(@(x) floor(-x.reward_start                 / cfg.timeBin), experiment.trial) + 1;
    trialStop         = arrayfun(@(x) floor((x.duration-1 - x.reward_start) / cfg.timeBin), experiment.trial) - 1;
    trialRange        = binOffset + (max(trialStart):min(trialStop));

%       longfigure; imagesc(alignedActivity(:,trialRange)); colorbar

    alignedActivity   = alignedActivity(:,trialRange,:);
    alignedTime       = alignedTime(trialRange);
    
    %%
    peakActivity      = nan(size(alignedActivity,1), 2, cfg.nShuffles + 1);
    for iShuffle = 0:cfg.nShuffles
      %% Use mean activity across trials to define time of maximum/minimum change from baseline
      baselineRate    = alignedActivity(:,alignedTime < 0,iShuffle+1);
      meanActivity    = mean(alignedActivity(:,:,iShuffle+1), 1);
      meanSigni       = ( meanActivity - mean(baselineRate(:)) ) / std(baselineRate(:));
      [~,maxTime]     = max(meanSigni);
      [~,minTime]     = min(meanSigni);
      
%       figure; plot(alignedTime, meanSigni)
      
      %% Distribution of deviations from baseline rates 
      peakActivity(:,:,iShuffle+1)  = ( alignedActivity(:,[maxTime, minTime],iShuffle+1) - mean(baselineRate,2)) ./ std(baselineRate,0,2);
%       peakActivity(:,:,iShuffle+1)  = alignedActivity(:,[maxTime, minTime],iShuffle+1) - mean(baselineRate,2);
    end
    
    %%
    peakActivePValue  = [ 1 + sum(peakActivity(:,1,1) > peakActivity(:,1,2:end), 3)     ...
                        , 1 + sum(peakActivity(:,2,1) < peakActivity(:,2,2:end), 3)     ...
                        ] ./ ( 1 + cfg.nShuffles );
    
    %% Unused : 2-tailed p-value for whether the activity in a given trial is more extreme that that of its shuffled versions
%     %{
    tailPValue        = 2 * min( ( 1 + sum(alignedActivity(:,:,1) < alignedActivity(:,:,2:end), 3) ) / ( 1 + cfg.nShuffles )   ...
                               , ( 1 + sum(alignedActivity(:,:,1) > alignedActivity(:,:,2:end), 3) ) / ( 1 + cfg.nShuffles )   ...
                               );
%     figure; imagesc(tailPValue < 0.05); colorbar
%     figure; bandplot(alignedTime,alignedActivity(:,:,1)'); hold on; bandplot(alignedTime, remould(permute(alignedActivity(:,:,2:end),[2 1 3]),2))
    
    %
    [isSignificant, pValueThreshold] = accumfun(1, @(x) fdrBenjaminiHochberg(tailPValue(x,:), 0.05), 1:numel(experiment.trial));
    figure; imagesc(isSignificant); colorbar
    figure; bandplot(alignedTime,alignedActivity(:,:,1)'); hold on; bandplot(alignedTime, remould(permute(alignedActivity(:,:,2:end),[2 1 3]),2))
    %}
  end
  
  %% Save models
  save(outputFile, 'decoder', 'cfg');
  fprintf(' -->  %s\n', outputFile);
  
end
