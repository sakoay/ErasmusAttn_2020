function plotFiringRatePrediction(modelFile)
  
  %% Load data
  if isstruct(modelFile)
    model             = modelFile;
  else
    model             = load(modelFile);
  end
  
  %% Diagnostic plots of model fits, per cell
  for iCell = 1:numel(model.categoryModel)
    if isempty(model.categoryModel{iCell})
      continue
    end
    
    %% Plot trial-average firing rate and predictions separately for all unique data conditions
    plotPredictionByCondition(model.categoryModel{iCell}, 'Fully specialized model', model.cfg);
    plotPredictionByCondition(model.hierarchicalModel{iCell}, 'Hierarchical model', model.cfg);
  end
  
end

function specs = parseSpecs(specs, lookup)
  specs         = regexp(specs, '([0-9]+)\s*=\s*([^,]+)', 'tokens');
  specs         = catcell(1, specs, 2);
  
  if isempty(specs)
    specs       = arrayfun(@num2str, lookup, 'uniformoutput', false);
  else
    specs       = specs{:};
    [~,iMatch]  = ismember(lookup, cellfun(@str2num,specs(:,1)));
    specs       = specs(iMatch,2);
  end
end

function fig = plotPredictionByCondition(model, description, cfg)
  
  %% Use the basic unspecialized model as a comparison
  refModel                = unique(model.firstAncestor);
  experiment              = refModel.design.dspec.expt;
  
  %% Specify categories of trials for which to plot data together
  trialConditions         = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);
  [trialCategory, ~, categIndex]    = unique(trialConditions, 'rows');
  
  %% Labels for trial categories
  arrayfun(@(x) parseSpecs(experiment.desc.(cfg.behavConditions{x}),trialCategory(:,x)), 1:numel(cfg.behavConditions), 'UniformOutput', false)
  
  %%
  CdirLabel               = parseSpecs(experiment.desc.direction_C   , config(:,1));
  condLabel               = parseSpecs(experiment.desc.condition_code, condition);
  
  %% Configure plots
  cellLabel               = [experiment.monkey ' ' experiment.cell_id];
  [pan,shape,fig]         = makePanels(size(trialCategory,1), experiment.id, experiment.id, 'aspectratio', 1, 'plotmargins', struct('b',20), 'panelmargins', struct('b',15,'t',15));
  
  %% Loop through categories of trials and plot each one
  for iCond = 1:numel(condition)
    for iConfig = 1:size(config,1)
      %% Select all trials within this category
      trial               = experiment.trials(condIndex == iCond & configIndex == iConfig);
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
                          = ksdensity(trial(iTrial).SS, time(trialRange), 'Bandwidth', 50);
        end
      end
      
      % Truncate data at last nonzero entry, for compactness of plotting
%       firingRate(find(sum(firingRate,2) > 0, 1, 'last') + 1:end,:)  = [];
%       time(size(firingRate,1) + 1:end)                              = [];
      
      %% Axes formatting
      axs                 = selectPanel(pan, [iCond,iConfig], shape);
      hold(axs, 'on');
      axis(axs, 'tight');
      colormap(axs, flipud(gray));
      if iCond == 1
        title(axs, sprintf('%s C, gap %.3g\\circ', CdirLabel{iConfig}, config(iConfig,2)), 'fontweight', 'normal');
      end
      if iConfig == 1
        ylabel(axs, condLabel{iCond});
      end
%       ylabel(axs, 'Firing rate (A.U.)');
      set(axs, 'xcolor', [1 1 1]*0.5, 'ycolor', [1 1 1]*0.5, 'clim', [0 1]);
      
      %% Plot firing rate with one trial per row
      rateRange           = quantile(firingRate(firingRate ~= 0), 0.95);
      yTrial              = (0:numel(trial) - 1) * rateRange;
      plot(axs, time, firingRate + yTrial, 'k');

      %% Indicate stimulus periods with a line
      yOffset             = numel(trial) * rateRange * 0.05;
      iOffset             = 0;
      ticks               = [];
      tickLabels          = {};
      for what = {'cue', 'c'}
        onTime            = zeros(numel(time), numel(trial));
        for iTrial = 1:numel(trial)
          if isnan(trial(iTrial).([what{:} '_start']))
            continue
          end
          onTime( 1 + (trial(iTrial).([what{:} '_start']):trial(iTrial).([what{:} '_stop'])), iTrial )  = 1;
        end
        onTime            = mean(onTime, 2);
        
        iOffset           = iOffset + 1;
        ticks(end+1)      = -iOffset * yOffset;
        tickLabels{end+1} = what{:};
        image('parent', axs, 'xdata', time([1 end]), 'ydata', -(iOffset + 0.15*[-1 1])*yOffset, 'cdata', onTime', 'CDataMapping', 'scaled');
      end
      set(axs, 'ytick', flip(ticks), 'yticklabel', flip(tickLabels));
      
      %% Indicate various behavioral events
      behavEvent          = {'fix_on', 'delay_start', 'choice_target', 'reward_start'};
      eventColor          = lines(numel(behavEvent));
      ticks               = [];
      tickLabels          = {};
      for iEvent = 1:numel(behavEvent)
        ticks(end+1)      = median([trial.(behavEvent{iEvent})]);
        tickLabels{end+1} = coloredText(regexprep(behavEvent{iEvent}, '_.+', ''), eventColor(iEvent,:));
        line('parent', axs, 'xdata', [trial.(behavEvent{iEvent})], 'ydata', yTrial, 'marker', '.', 'markersize', 10, 'linestyle', 'none', 'color', eventColor(iEvent,:));
      end
      
      if any(isnan(ticks))
        xlabel(axs, 'Time (ms)');
      else
        [ticks,iOrder]    = sort(ticks);
        set(axs, 'xtick', ticks, 'xticklabel', tickLabels(iOrder), 'xticklabelrotation', 90);
      end
      
      %% Use grid to indicate time scale
      drawnow;
      xRange              = get(axs, 'xlim');
      yRange              = get(axs, 'ylim');
      xGrid               = 0:1000:xRange(2);         % 1 second intervals
      uistack(line( 'parent', axs, 'xdata', colvec(repmat(xGrid,3,1)), 'ydata', colvec(repmat([yRange(:); nan],1,numel(xGrid)))       ...
                  , 'linestyle', '-.', 'color', [1 1 1]*0.7 ), 'bottom');
    end
  end
  
end
