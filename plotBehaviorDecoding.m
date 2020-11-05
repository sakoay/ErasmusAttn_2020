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

  %% Plots for each decoded quantity
%   %{
  [pan, shape, fig]   = makePanels( [numel(cfg.behavCategories),numel(cfg.behavEvents)], cfg.animal, cfg.animal     ...
                                  , 'aspectratio', 1.25, 'twosided', true, 'maxsubplotsize', 200 );
  for iBehav = 1:numel(cfg.behavCategories)
    plotDecodingStatistics(pan, iBehav, shape, cat(1,decoder.(cfg.behavCategories{iBehav})), 2, cfg);
  end
  %}
  
  %% Plots for decoders that are further separated by condition values
%   %{
  nSubsets            = cellfun(@(x) numel(decoder(1).(x)), cfg.behavCategories);
  [pan, shape, fig]   = makePanels( [sum(nSubsets(nSubsets > 1)),numel(cfg.behavEvents)], cfg.animal, cfg.animal     ...
                                  , 'aspectratio', 1.25, 'twosided', true, 'maxsubplotsize', 200 );
  iDecode             = 1;
  for iBehav = 1:numel(cfg.behavCategories)
    if nSubsets(iBehav) > 1
      iDecode         = plotDecodingStatistics(pan, iDecode, shape, cat(1,decoder.(cfg.behavCategories{iBehav})), false, cfg);
    end
  end
  %}
  
end

function iRow = plotDecodingStatistics(pan, startRow, shape, decoder, averageAccuracy, cfg)
  
  %% Decoding using neural data aligned to different behavioral events
  for iAlign = 1:numel(cfg.behavEvents)
    %% Collect results along a common time axis across cells
    timeRange         = accumfun(2, @(x) extrema(x.model(iAlign).alignedTime(:)), decoder);
    timeRange         = [min(timeRange(1,:)), max(timeRange(2,:))];
    alignedTime       = timeRange(1):cfg.timeBin:timeRange(2);
    
    accuracy          = nan([numel(alignedTime), numel(decoder), cfg.nShuffles+1]);
    for iModel = 1:numel(decoder)
      model           = decoder(iModel).model(iAlign);
      iTime           = binarySearch(alignedTime, model.alignedTime, 0, 0);
      accuracy(iTime,iModel,:)  = model.accuracy;
    end
    accuracy          = reshape(accuracy, [numel(alignedTime), size(decoder), cfg.nShuffles+1]);
    
    %% Optionally compute average accuracy across decoders for different subsets of condition values
    switch averageAccuracy
      case 0
        %% No averaging
        pValue        = mean(accuracy(:,:,:,2:end) >= accuracy(:,:,:,1), 4);
      case 1
        %% Average accuracies across decoders, and then compute p-value per time bin and neuron
        accuracy      = mean(accuracy, 3);
        pValue        = mean(accuracy(:,:,:,2:end) >= accuracy(:,:,:,1), 4);
      case 2
        %% Use best p-value across subset-specific decoders
        pValue        = mean(accuracy(:,:,:,2:end) >= accuracy(:,:,:,1), 4);
        pValue(isnan(accuracy(:,:,:,1)))  = nan;
        nModels       = sum(~isnan(pValue), 3);
        pValue        = nModels .* min(pValue,[],3);    % Bonferroni correction for taking the best of two tests
        
        %% Display average accuracy across decoders
        accuracy      = mean(accuracy, 3);
      otherwise
        error('plotDecodingStatistics:averageAccuracy', 'averageAccuracy must be either true, false or 2 for computing best p-value across decoder types.');
    end
    nullAccuracy      = mean(accuracy(:,:,:,2:end), 4);
    accuracy(:,:,:,2:end)   = [];
    pValue(isnan(accuracy)) = nan;
    
    %% Plot separately for each decoder
    iRow              = startRow;
    for iDecode = 1:size(accuracy,3)
      %% Define threshold for statistical significance after correcting for false discovery rate
      [isSignificant, pValueThreshold]      ...
                      = fdrBenjaminiHochberg(pValue(:,:,iDecode), cfg.significance);
%       isSignificant   = pValue(:,:,iDecode) <= cfg.significance;
      [phat, pci]     = binointerval(sum(isSignificant,2), sum(~isnan(pValue(:,:,iDecode)),2), cfg.significance);

      %% Decoding accuracy across cells
      axs             = selectPanel(pan, [iRow,iAlign], shape);
      yRange          = [0, 1];
%       yRange          = [0, max(accuracy(:,:,iDecode),[],'all')];
      xlabel(axs, sprintf('Time from %s', strrep(model.alignment,'_',' ')));
      ylabel(axs, sprintf('%s accuracy', strrep(decoder(1).variable,'_',' ')));
      set(axs, 'xlim', alignedTime([1 end]), 'ylim', yRange);
      if averageAccuracy
        title(axs, sprintf('%s = %s', strrep(decoder(1,iDecode).variable,'_',' '), mat2str(decoder(1,iDecode).values)), 'FontWeight', 'normal');
      end

      line('parent', axs, 'xdata', [0 0], 'ydata', yRange, 'linewidth', 1.5, 'linestyle', ':', 'color', [1 1 1]*0.7);
      line('parent', axs, 'xdata', alignedTime, 'ydata', nanmean(nullAccuracy(:,:,iDecode),2), 'linewidth', 1.5, 'linestyle', '-.', 'color', [1 1 1]*0.5);
      bandplot(axs, alignedTime, accuracy(:,:,iDecode), cfg.confIntFcn);
%     text(alignedTime(end), nanmean(nullAccuracy(:)), ' (chance)', 'parent', axs, 'fontsize', 14, 'color', [1 1 1]*0.5, 'clipping', 'off');

      %% Fraction of cells with significant decoding
      yyaxis(axs, 'right');
      ylabel(axs, '% cells significant', 'rotation', -90, 'verticalalignment', 'bottom');
      set(axs, 'ylim', [0 100]);
      color           = get(axs, 'ycolor');

      patch('parent', axs, 'xdata', [alignedTime(:); flip(alignedTime(:))], 'ydata', 100*[pci(:,1); flip(pci(:,2))], 'facecolor', color, 'facealpha', 0.25, 'linestyle', 'none');
      line('parent', axs, 'xdata', alignedTime, 'ydata', 100*phat, 'linewidth', 1.5, 'linestyle', '-', 'color', color, 'clipping', 'off');
      
      %% Increment row for plots
      iRow            = iRow + 1;
    end
  end
  
end

