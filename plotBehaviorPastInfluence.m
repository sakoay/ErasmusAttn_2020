function plotBehaviorPastInfluence(modelFile)
  
  %% Load data
  if ischar(modelFile)
    model             = load(modelFile);
  else
    model             = modelFile;
  end
  
  %% Variable names
  Xvariables          = model.behaviorModel{1}.Xvariables;
  isPastVar           = strncmp(Xvariables, 'past', 4);
  Xvariables          = regexprep(Xvariables, 'saccade', 'sacc');
  Xvariables          = regexprep(Xvariables, 'gap(.+)', 'gap^{[$1]}');
  Xvariables          = regexprep(Xvariables, '^past(.+)_(.+)', '$2_{(-$1)}');
  
  taskLabel           = {'L/R task', 'U/D task'};
  
  %% Compare fits across animals
  aniIndex            = find(~contains(model.experiment, '_')');
  modelColor          = lines(size(model.behaviorModel,2));
  
  [pan,shape,fig]     = makePanels([1,2*numel(aniIndex)], '', '', 'aspectratio', 0.6, 'maxsubplotsize', 250, 'plotmargins', struct('l',30), 'panelmargins', struct('b',25));
  for iAni = 1:numel(aniIndex)
    for iPast = 0:1
      selVar          = isPastVar == iPast;
      varIndex        = 1:sum(selVar);
      xRange          = varIndex([1 end]) + [-0.5,0.5];
      
      %% Axes formatting
      axs             = selectPanel(pan, 1 + 2*(iAni-1) + iPast, shape);
      hold(axs, 'on');
      set(axs, 'xlim', xRange, 'xtick', varIndex, 'xticklabel', Xvariables(selVar), 'xticklabelrotation', 90, 'fontsize', 16, 'fontsizemode', 'manual');
      
      if iPast
        title(axs, model.experiment{aniIndex(iAni)}, 'FontWeight', 'normal');
        nudge(axs, [-0.15 0 0.15 0]);
      else
        nudge(axs, [0.01 0 -0.1 0]);
        ylabel(axs, 'Weight to predict saccade');
%         set(axs, 'ylim', [0 2]);
      end
      line('parent', axs, 'xdata', xRange, 'ydata', [0 0], 'linestyle', '-.', 'linewidth', 1, 'color', [1 1 1]*0.7);

      %% Regressor weights for horizontal vs. vertical categories of trials
      hAni            = gobjects(size(aniIndex));
      for iModel = 1:size(model.behaviorModel,2)
        %% Model weights and confidence intervals for pooled data
        bootstrapW    = model.behaviorModel{aniIndex(iAni),iModel}.bootstrapW(selVar,:);
        modelW        = mean(bootstrapW, 2);
        ciW           = quantile(bootstrapW - modelW, [0.05,0.95], 2);
        hAni(iModel)  = errorbar(axs, varIndex - 0.05*(2*iModel-3), modelW, -ciW(:,1), ciW(:,2), 'linestyle', 'none', 'color', modelColor(iModel,:), 'marker', 'o', 'markersize', 5, 'linewidth', 1.5);
      end
    end
    
    legend(hAni, taskLabel, 'location', 'best', 'box', 'on');
  end
  
  %% Separately plot fits for each animal
  for iAni = 1:numel(aniIndex)
    %% Model with data per cell for a given animal
    pooledModel       = model.behaviorModel(aniIndex(iAni),:);
    selCell           = strncmp(model.experiment, [model.experiment{aniIndex(iAni)} '_'], numel(model.experiment{aniIndex(iAni)})+1);
    cellModel         = model.behaviorModel(selCell,:);
    
    %% Unused : Plots for horizontal vs. vertical trial models
    %{
    [pan,shape,fig]   = makePanels([1,size(model.behaviorModel,2)], '', '', 'aspectratio', 1.2, 'maxsubplotsize', 300, 'panelmargins', struct('b',25));
    for iModel = 1:size(model.behaviorModel,2)
      varIndex        = 1:sum(isPastVar);
      xRange          = varIndex([1 end]) + [-0.5,0.5];

      %% Axes formatting
      axs             = selectPanel(pan, iModel, shape);
      hold(axs, 'on');
      set(axs, 'xlim', xRange, 'xtick', varIndex, 'xticklabel', Xvariables(selVar), 'xticklabelrotation', 90, 'fontsize', 16, 'fontsizemode', 'manual');
      title(axs, [model.experiment{aniIndex(iAni)} ' : ' taskDir{iModel} ' task'], 'FontWeight', 'normal');
      line('parent', axs, 'xdata', xRange, 'ydata', [0 0], 'linestyle', '-.', 'linewidth', 1, 'color', [1 1 1]*0.7);

      %% Model weights for data in individual cell recording sessions
      selModel        = cellModel(:,iModel);
      selModel        = selModel(~cellfun(@isempty,selModel));
      bootstrapW      = accumfun(2, @(x) mean(x.bootstrapW( selVar,:),2), selModel);
      refW            = accumfun(2, @(x) mean(x.bootstrapW(~selVar,:),2), selModel);
      line( 'parent', axs, 'xdata', colvec( repmat(varIndex(:) - 0.05*(2*iModel-3),1,size(bootstrapW,2)) + 0.1*randn(size(bootstrapW)) )   ...
          , 'ydata', colvec(bootstrapW ./ refW(iModel,:)), 'linestyle', 'none', 'color', [1 1 1]*0.5, 'marker', '.', 'markersize', 5, 'linestyle', 'none' );

      %% Model weights and confidence intervals for pooled data
      refW            = pooledModel{iModel}.bootstrapW(~selVar,:);
      bootstrapW      = pooledModel{iModel}.bootstrapW( selVar,:) ./ refW(iModel,:);
      modelW          = mean(bootstrapW ./ refW(iModel,:), 2);
      ciW             = quantile(bootstrapW - modelW, [0.05,0.95], 2);
      hAni            = errorbar(axs, varIndex - 0.05*(2*iModel-3), modelW, -ciW(:,1), ciW(:,2), 'linestyle', 'none', 'color', [0 0 0], 'marker', 'o', 'markersize', 5, 'linewidth', 1.5);
    end
    %}

    %% Plots for horizontal vs. vertical trial models
%     %{
    [pan,shape,fig]   = makePanels([1,2*size(model.behaviorModel,2)], '', '', 'aspectratio', 0.5, 'maxsubplotsize', 300, 'plotmargins', struct('l',30), 'panelmargins', struct('b',25));
    for iModel = 1:size(model.behaviorModel,2)
      for iPast = 0:1
        selVar        = isPastVar == iPast;
        varIndex      = 1:sum(selVar);
        xRange        = varIndex([1 end]) + [-0.5,0.5];

        %% Axes formatting
        axs           = selectPanel(pan, 1 + 2*(iModel-1) + iPast, shape);
        hold(axs, 'on');
        axis(axs, 'tight');
        set(axs, 'xlim', xRange, 'xtick', varIndex, 'xticklabel', Xvariables(selVar), 'xticklabelrotation', 90, 'fontsize', 16, 'fontsizemode', 'manual');

        if iPast
          title(axs, [model.experiment{aniIndex(iAni)} ' : ' taskLabel{iModel}], 'FontWeight', 'normal');
          nudge(axs, [-0.15 0 0.15 0]);
        else
          nudge(axs, [0.01 0 -0.1 0]);
          ylabel(axs, 'Weight to predict saccade');
        end
        line('parent', axs, 'xdata', xRange, 'ydata', [0 0], 'linestyle', '-.', 'linewidth', 1, 'color', [1 1 1]*0.7);

        %% Model weights for data in individual cell recording sessions
        selModel      = cellModel(:,iModel);
        bootstrapW    = accumfun(2, @(x) mean(x.bootstrapW(selVar,:),2), selModel(~cellfun(@isempty,selModel)));
        line( 'parent', axs, 'xdata', colvec( repmat(varIndex(:) - 0.05*(2*iModel-3),1,size(bootstrapW,2)) + 0.1*randn(size(bootstrapW)) )   ...
            , 'ydata', bootstrapW(:), 'linestyle', 'none', 'color', [1 1 1]*0.5, 'marker', '.', 'markersize', 5, 'linestyle', 'none' );
        
        %% Model weights and confidence intervals for pooled data
        bootstrapW    = pooledModel{iModel}.bootstrapW(selVar,:);
        modelW        = mean(bootstrapW, 2);
        ciW           = quantile(bootstrapW - modelW, [0.05,0.95], 2);
        hAni          = errorbar(axs, varIndex - 0.05*(2*iModel-3), modelW, -ciW(:,1), ciW(:,2), 'linestyle', 'none', 'color', [0 0 0], 'marker', 'o', 'markersize', 5, 'linewidth', 1.5);
      end
    end
    %}
  end
  
end
