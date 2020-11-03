function plotBehaviorDecoding(modelFile)
  
  %% Load data
  load(modelFile);
  
  % Omit cells that could not be fit due to data restrictions
  decoder(~cellfun(@(x) isfield(x,'experiment'), decoder)) = [];
  decoder             = [decoder{:}];
  
  %% Analysis configuration
  cfg.animal          = regexprep(decoder(1).experiment.id, '_.*', '');

  cfg.confIntFcn      = [normcdf([-1 1],0,1); 0.025 0.975];
  cfg.significance    = 0.05;

  %% Configure plots
  [pan, shape, fig]   = makePanels( [numel(cfg.behavCategories),numel(cfg.behavEvents)], cfg.animal, cfg.animal     ...
                                  , 'aspectratio', 1.5, 'twosided', true );
  
  %% Plots for each decoded quantity
  for iBehav = 1:numel(cfg.behavCategories)
    plotDecodingStatistics(pan, iBehav, shape, [decoder.(cfg.behavCategories{iBehav})], cfg.behavCategories{iBehav}, cfg);
  end
  
end

function plotDecodingStatistics(pan, iRow, shape, decoder, description, cfg)
  
  %% Decoding using neural data aligned to different behavioral events
  for iAlign = 1:numel(cfg.behavEvents)
    %% Collect results along a common time axis across cells
    timeRange         = accumfun(2, @(x) extrema(x.model(iAlign).alignedTime(:)), decoder);
    timeRange         = [min(timeRange(1,:)), max(timeRange(2,:))];
    alignedTime       = timeRange(1):cfg.timeBin:timeRange(2);
    
    accuracy          = nan(numel(alignedTime), numel(decoder));
    nullAccuracy      = nan(numel(alignedTime), numel(decoder));
    pValue            = nan(numel(alignedTime), numel(decoder));
    for iModel = 1:numel(decoder)
      model           = decoder(iModel).model(iAlign);
      iTime           = binarySearch(alignedTime, model.alignedTime, 0, 0);
      accuracy(iTime,iModel)      = model.accuracy(:,1);
      nullAccuracy(iTime,iModel)  = mean(model.accuracy(:,2:end), 2);
      pValue(iTime,iModel)        = mean(model.accuracy(:,2:end) >= model.accuracy(:,1), 2);     % N.B. a pseudocount is conservatively included for the unshuffled experiment
    end
    
    %% Define threshold for statistical significance after correcting for false discovery rate
    [isSignificant, pValueThreshold]      ...
                      = fdrBenjaminiHochberg(pValue, cfg.significance);
%     isSignificant     = pValue <= cfg.significance;
    
    %% Decoding accuracy across cells
    axs               = selectPanel(pan, [iRow,iAlign], shape);
    yRange            = [0, max(accuracy(:))];
    xlabel(axs, sprintf('Time from %s', strrep(model.alignment,'_',' ')));
    ylabel(axs, sprintf('%s accuracy', strrep(description,'_',' ')));
    set(axs, 'xlim', alignedTime([1 end]), 'ylim', yRange);

    line('parent', axs, 'xdata', [0 0], 'ydata', yRange, 'linewidth', 1.5, 'linestyle', ':', 'color', [1 1 1]*0.7);
    line('parent', axs, 'xdata', alignedTime, 'ydata', nanmean(nullAccuracy,2), 'linewidth', 1.5, 'linestyle', '-.', 'color', [1 1 1]*0.5);
    bandplot(axs, alignedTime, accuracy, cfg.confIntFcn);
%     text(alignedTime(end), nanmean(nullAccuracy(:)), ' (chance)', 'parent', axs, 'fontsize', 14, 'color', [1 1 1]*0.5, 'clipping', 'off');

    %% Fraction of cells with significant decoding
    [phat, pci]       = binointerval(sum(isSignificant,2), sum(~isnan(pValue),2), cfg.significance);
    yyaxis(axs, 'right');
    ylabel(axs, '% cells significant', 'rotation', -90, 'verticalalignment', 'bottom');
    set(axs, 'ylim', [0 100]);
    color             = get(axs, 'ycolor');
    
    patch('parent', axs, 'xdata', [alignedTime(:); flip(alignedTime(:))], 'ydata', 100*[pci(:,1); flip(pci(:,2))], 'facecolor', color, 'facealpha', 0.25, 'linestyle', 'none');
    line('parent', axs, 'xdata', alignedTime, 'ydata', 100*phat, 'linewidth', 1.5, 'linestyle', '-', 'color', color, 'clipping', 'off');
  end
  
end

