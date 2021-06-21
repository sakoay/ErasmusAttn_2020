function plotCorrectVsErrorTrials(dataFile)
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
  end
  
  %% Plot some quasi-random sample of cells
%   for iCell = 1:4:numel(data.cellData)
  for iCell = 1:numel(data.cellData)
    plotCellData(data.cellData(iCell), data.experiment);
  end
  
end

function fig = plotCellData(data, experiment)
  
  %% Select trials of interest
%   data.trials(~ismember([data.trials.past_condition], 1:2)) = [];
  data.trials(~ismember([data.trials.condition_code], 4)) = [];
  if isempty(data.trials)
    return;
  end
  
  %% Plotting configuration
  gapSignMap              = [  270     -1     ... left
                            ;   90     +1     ... right
                            ;  180     -1     ... down
                            ;    0     +1     ... up
                            ];
  gapDirLabel             = {'L/R', 'U/D'};
                      
  %% Specify categories of trials for which to plot data together
  [~,trialConfig]         = ismember(cat(1,data.trials.gap_direction), gapSignMap(:,1));
  [config,~,configIndex]  = unique(ceil(trialConfig/2), 'rows');
  [condition,~,condIndex] = unique([data.trials.condition_code]);
  
  %% Labels for trial categories
  condLabel               = parseVariableSpecifications(experiment.desc.condition_code, config(:,1), true);
  configLabel             = gapDirLabel(config);
  
  %% Past info categories
  [pastConfig,~,pastIndex]= unique([cat(1,data.trials.past_gap), cat(1,data.trials.past_saccade)], 'rows');
  
  %% Configure plots
  cellLabel               = [data.monkey ' ' data.cell_id];
  [pan,shape,fig]         = makePanels([numel(condition),size(config,1)], cellLabel, cellLabel, 'aspectratio', 1, 'plotmargins', struct('b',20), 'panelmargins', struct('b',15,'t',15));
  
  %% Loop through categories of trials and plot each one
  for iCond = 1:numel(condition)
    for iConfig = 1:size(config,1)
      %% Select all trials within this category
      selTrials           = condIndex == iCond & configIndex == iConfig;
      trial               = data.trials(selTrials);
      if isempty(trial)
        continue
      end
      
      %% Convert spike trains to smoothed firing rates
      time                = 0:max([trial.duration]);
      firingRate          = nan(numel(time), numel(trial));
      for iTrial = 1:numel(trial)
        trialRange        = 1:trial(iTrial).duration;
        if isempty(trial(iTrial).SS)
          firingRate(trialRange,iTrial)	= 0;
        else
          firingRate(trialRange,iTrial)   ...
                          = ksdensity(trial(iTrial).SS, time(trialRange), 'Bandwidth', 100);
        end
      end
      
      % Truncate data at last nonzero entry, for compactness of plotting
%       firingRate(find(sum(firingRate,2) > 0, 1, 'last') + 1:end,:)  = [];
%       time(size(firingRate,1) + 1:end)                              = [];
      
      %% Axes formatting
      axs                 = selectPanel(pan, [iCond,iConfig], shape);
      hold(axs, 'on');
      axis(axs, 'tight');
      title(axs, [configLabel{iConfig} ' ' condLabel{iCond}], 'fontweight', 'normal');
      set(axs, 'xcolor', [1 1 1]*0.5, 'ycolor', [1 1 1]*0.5, 'clim', [0 1]);
      
      %% Plot firing rate for various past trial configurations
      selPast             = pastIndex(selTrials);
      for iPast = 1:size(pastConfig,1)
        bandplot(axs, time, firingRate(:,selPast == iPast));
      end
    end
  end
  
end
