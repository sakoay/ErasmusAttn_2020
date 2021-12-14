function erasmus_decoding_figure(modelFile)

  %% Figure setup
%   PaneledFigure([4 4], 'small')
  figMain       = PaneledFigure(reshape(1:20,5,[])', [160 190], [15 5 20 0], [2 -30]);
  figMain.newGroup((2:5) + (0:5:15)') = false;
  iMain         = 0;
  
  %% Load data
  model                     = cell(size(modelFile));
  for iFile = 1:numel(modelFile)
    model{iFile}            = load(modelFile(iFile).name);
    cfg                     = model{iFile}.cfg;
    model{iFile}            = model{iFile}.decoder;
    
    %% Omit cells that could not be fit due to data restrictions
    model{iFile}(~cellfun(@(x) isfield(x,'experiment'), model{iFile}))  = [];
    model{iFile}            = [model{iFile}{:}];
    
    %% Ensure that all decoders have the same fields
    for iBehav = 1:numel(cfg.behavCategories)
      for iDec = 1:numel(model{iFile})
        assert(isfield(model{iFile}(iDec).(cfg.behavCategories{iBehav}), 'model'));
      end
    end
  end
  
  %% Decoding accuracy for various event alignments
  iMain         = plotDecodability(figMain, iMain, model, cfg);
  
  %% Save figures
  codeFile      = mfilename('fullpath');
  figMain.export(codeFile, [], [], [], repmat(0.12 + 0.18j, figMain.size));
  
end

function iPlot = plotDecodability(fig, iPlot, aniModel, cfg)

  %% Analysis configuration
  aniColor                = [FormatDefaults.darkBlue; FormatDefaults.mediumPurple];
  cfg.confIntFcn          = [normcdf([-1 1],0,1); 0.025 0.975];
  cfg.significance        = 0.05;
  xRange                  = [ 0     3     ... fix_on
                            ; nan   nan   ... cue_start
                            ; -0.5  2.5   ... c_start
                            ; -1    2     ... delay_start
                            ; nan   nan   ... choice_target
                            ; -1.5  1.5   ... saccade_onset
                            ; -2    1     ... reward_start
                            ];
  
  %% Decoding using neural data aligned to different behavioral events
  for iBehav = 1:numel(cfg.behavCategories)
    iStart                = iPlot;
    hPlot                 = gobjects(numel(aniModel), numel(cfg.behavEvents));
    for iAni = 1:numel(aniModel)
      %% Collect decoders of a given type across all cells recorded for this animal
      decoder             = cat(1, aniModel{iAni}.(cfg.behavCategories{iBehav}));
      decodeLabel         = strrep(strrep(regexprep(strrep(strrep(cfg.behavCategories{iBehav},'_',' '), 'condition', 'outcome'), ' (code|direction)', ''), '-condition', '-cond'), 'past ', 'past-');
      iPlot               = iStart;

      %% Plot decoding measures aligned to various behavioral events
      for iAlign = 1:numel(cfg.behavEvents)
        if ~isempty(regexp(cfg.behavEvents{iAlign}, 'cue|choice', 'match', 'once'))
          continue
        end
        
        %% Collect results along a common time axis across cells
        timeRange         = accumfun(2, @(x) extrema(x.model(iAlign).alignedTime(:)), decoder(~arrayfun(@(x) isempty(x.model), decoder)));
        timeRange         = [min(timeRange(1,:)), max(timeRange(2,:))];
        alignedTime       = timeRange(1):cfg.timeBin:timeRange(2);

        accuracy          = nan([numel(alignedTime), numel(decoder), cfg.nShuffles+1]);
        for iModel = 1:numel(decoder)
          if isempty(decoder(iModel).model)
            continue
          end
          model           = decoder(iModel).model(iAlign);
          iTime           = binarySearch(alignedTime, model.alignedTime, 0, 0);
          accuracy(iTime,iModel,:)        = model.accuracy;
        end
        accuracy          = reshape(accuracy, [numel(alignedTime), size(decoder), cfg.nShuffles+1]);

        %% Use best p-value across decoders for different subsets of condition values
        pValue            = mean(accuracy(:,:,:,2:end) >= accuracy(:,:,:,1), 4);
        pValue(isnan(accuracy(:,:,:,1)))  = nan;
        nModels           = sum(~isnan(pValue), 3);
        pValue            = nModels .* min(pValue,[],3);    % Bonferroni correction for taking the best of two tests
        
        accuracy(:,:,:,2:end)             = [];
        pValue(isnan(all(accuracy,2)))    = nan;

        %% Define threshold for statistical significance after correcting for false discovery rate
        [isSignificant, pValueThreshold]      ...
                          = fdrBenjaminiHochberg(pValue, cfg.significance);
        [phat, pci]       = binointerval(sum(isSignificant,2), sum(~isnan(pValue),2), cfg.significance);

        %% Fraction of cells with significant decoding
        [axs, iPlot]      = fig.panel(iPlot);
        xlabel(axs, sprintf('{\\itt} from %s (s)', regexprep(regexprep(strrep(strrep(model.alignment,'_',' '), 'fix on', 'fixation'), ' (start|onset)$', ''), '^c$', 'C onset')));
        if iAlign == 1
          if numel(decodeLabel) > 10
            tag           = 'dec.';
          else
            tag           = 'decoding';
          end
          ylabel(axs, multitext([], '% cells w/ signi.', [decodeLabel, ' ', tag]));
        else
          set(axs, 'ytick', [0 60]);
        end
        set(axs, 'ylim', [0 60], 'xlim', xRange(iAlign,:), 'clipping', 'on');

        patch('parent', axs, 'xdata', [alignedTime(:); flip(alignedTime(:))]/1000, 'ydata', 100*[pci(:,1); flip(pci(:,2))], 'facecolor', aniColor(iAni,:), 'facealpha', 0.25, 'linestyle', 'none');
        hPlot(iAni,iAlign)= line('parent', axs, 'xdata', alignedTime/1000, 'ydata', 100*phat, 'linewidth', 1.5, 'linestyle', '-', 'color', aniColor(iAni,:));
      end
    end
        
    %% Legend
    hLegend               = legend(hPlot(:,1), cellfun(@(x) regexprep(x(1).experiment.id, '_.*', ''), aniModel, 'UniformOutput', false), 'location', 'best', 'box', 'off');
  end
  
end
