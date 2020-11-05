function model = plotBehaviorDependence(modelFile)
  
  %% Load data
  if ischar(modelFile)
    modelFile         = {modelFile};
  end
  model               = cellfun(@load, modelFile, 'UniformOutput', false);

  %% Compare fits across animals
  aniModel            = accumfun(1, @(x) x.behaviorModel{end}, model);
  aniColor            = lines(numel(aniModel));
  varIndex            = 1:numel(aniModel(1).Xvariables);
  
  [pan,shape,fig]     = makePanels([1,numel(model)+1], '', '', 'aspectratio', 1.8, 'maxsubplotsize', 600, 'panelmargins', struct('b',40));
  axs                 = selectPanel(pan, 1, shape);
  hold(axs, 'on');
  set(axs, 'xlim', 0.5+[0, numel(aniModel(1).Xvariables)], 'xtick', varIndex, 'xticklabel', aniModel(1).Xvariables, 'xticklabelrotation', 90);
  ylabel(axs, 'Weight to predict trial outcome');
  line('parent', axs, 'xdata', 0.5+[0, numel(aniModel(1).Xvariables)], 'ydata', [0 0], 'linestyle', '-.', 'linewidth', 1, 'color', [1 1 1]*0.7);
  
  hAni                = gobjects(size(aniModel));
  for iAni = 1:numel(aniModel)
    %% Model weights and confidence intervals for pooled data
    fitW              = aniModel(iAni).modelW(:,aniModel(iAni).IndexMinDeviance);
    ciW               = quantile(aniModel(iAni).bootstrapW - fitW, [0.05,0.95], 2);
    hAni(iAni)        = errorbar(axs, varIndex - 0.05*(2*iAni-3), fitW, -ciW(:,1), ciW(:,2), 'linestyle', 'none', 'color', aniColor(iAni,:), 'marker', 'o', 'markersize', 5, 'linewidth', 1.5);
  end
  
  legend(hAni, cellfun(@(x) x.experiment{end}, model, 'UniformOutput', false), 'location', 'best');
  
  
  %% Separately plot fits for each animal
  for iModel = 1:numel(model)
    %% Model with data pooled from all cells 
    pooledModel       = model{iModel}.behaviorModel{end};
    
    %% Collect fit results across cells
    cellModel         = model{iModel}.behaviorModel(1:end-1);
    cellModel(cellfun(@isempty, cellModel)) = [];
    modelW            = nan(size(pooledModel.modelW,1), numel(cellModel));
    for iCell = 1:numel(cellModel)
      modelW(ismember(pooledModel.Xvariables,cellModel{iCell}.Xvariables),iCell)    ...
                      = cellModel{iCell}.modelW(:,cellModel{iCell}.IndexMinDeviance);
    end
    
    %% Axes and labels 
    axs               = selectPanel(pan, iModel+1, shape);
    hold(axs, 'on');
    set(axs, 'xlim', 0.5+[0, numel(pooledModel.Xvariables)], 'xtick', varIndex, 'xticklabel', pooledModel.Xvariables, 'xticklabelrotation', 90);
    ylabel(axs, 'Weight to predict trial outcome');
    title(axs, model{iModel}.experiment{end});
    line('parent', axs, 'xdata', 0.5+[0, numel(pooledModel.Xvariables)], 'ydata', [0 0], 'linestyle', '-.', 'linewidth', 1, 'color', [1 1 1]*0.7);
    
    %% Model weights for individual cells
    line('parent', axs, 'xdata', colvec(varIndex(:) + 0.1*randn(1,numel(cellModel))), 'ydata', modelW(:), 'marker', '.', 'markersize', 10, 'color', [1 1 1]*0.7, 'linestyle', 'none');
    
    %% Model weights and confidence intervals for pooled data
    fitW              = pooledModel.modelW(:,pooledModel.IndexMinDeviance);
    ciW               = quantile(pooledModel.bootstrapW - fitW, [0.05,0.95], 2);
    errorbar(axs, varIndex, fitW, -ciW(:,1), ciW(:,2), 'linestyle', 'none', 'color', [0.8 0 0], 'marker', 'o', 'markersize', 5, 'linewidth', 1.5);
  end
  
end
