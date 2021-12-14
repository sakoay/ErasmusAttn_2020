function erasmus_encoding_figure(modelFile)

  %% Figure setup
%   PaneledFigure([4 4], 'small')
  layout        = [   1   1   7   9  11  13     ...
                  ;   2   2   7  10  12  14     ...
                  ;   3   4   8  15  17  19     ...
                  ;   5   6   8  16  18  20     ...
                  ];
  labelOffset   = [ 0.35j   0   0     0.4-0.45j   0   0     ...
                  ;     0   0   0.12  0           0   0     ...
                  ;   0.3   0   0     0.4-0.45j   0   0     ...
                  ;     0   0   0.2   0           0   0     ...
                  ];
  figMain       = PaneledFigure(layout, [125 125], [10 10 10 10], [0 15]);
  figMain.newGroup([2 4:6 10:14 16:20]) = false;
  iMain         = 0;

  %% Custom panel sizes
  nudge(figMain.panel(1:2), [0 0 -0.6 0]);
  nudge(figMain.panel(1), [0 0.05]);
  nudge(figMain.panel(2), [-0.6 -0.1 0.6 0]);
  nudge(figMain.panel([3 5]), [-0.05 0]);
  nudge(figMain.panel([4 6]), [-0.6 0]);
  nudge(figMain.panel(7:8), [-0.7 0 1.2 -0.5]);
  nudge(figMain.panel(7), [0 0.45]);
  nudge(figMain.panel(8), [0 -0.5]);
  nudge(figMain.panel(9:20), [-0.15 0.05 0.3 0]);
  nudge(figMain.panel(9:14), [0 0.15]);
  for iRow = [0 6]
    nudge(figMain.panel(iRow + (9:2:13)), [0 -0.55]);
    nudge(figMain.panel(iRow + (9:10)), [0.6 0]);
  end
  for iRow = [0:1, 6:7]
    figMain.distributeHorizontally(iRow + (9:2:13));
  end
  
  %% Load data
  model                     = cell(size(modelFile));
  for iFile = 1:numel(modelFile)
    model{iFile}            = load(modelFile(iFile).name);
    
    %% Postprocessing 
    cfg                     = model{iFile}.cfg;
    catLabel                = regexprep(regexprep(cfg.behavConditions, '^(past)_(.)', '$1${upper($2)}'), '_.*', '');
    unspecializedModel      = model{iFile}.unspecializedModel;
    hierarchicalModel       = model{iFile}.hierarchicalModel;
    prediction              = cell(size(unspecializedModel));
    parfor iCell = 1:numel(unspecializedModel)
      prediction{iCell}     = modelTrialPredictions(unspecializedModel{iCell}, hierarchicalModel{iCell}, cfg, catLabel);
    end
    model{iFile}.prediction = [prediction{:}];
  end
  
  %% Model illustration
  iMain         = plotModelIlustration(figMain, iMain, model{1}.unspecializedModel{1});
  
  %% Modulation strengths
  iMain         = plotModulationStrength(figMain, iMain, model);
  
  %% Specialization categories
  iMain         = plotModulationCategories(figMain, iMain, model);
  
  %% Model predictions for example cells
  iMain         = plotModelPredictions(figMain, iMain, model);
  
  %% Save figures
  codeFile                = mfilename('fullpath');
  figMain.export(codeFile, [], [], [], labelOffset);
  
end

function combo = modelTrialPredictions(refModel, model, cfg, categoryLabels)

  %% Get frames corresponding to each trial
  assert(numel(refModel) == 1);

  [trialID, frames]   = SplitVec(refModel.design.trialX, 'equal', 'firstval', 'index');
%   trialEpochs       = cellfun(@(x) refModel.design.dspec.expt.eventEpochs(x), frames, 'UniformOutput', false);

  %% Preallocate output
  combo               = mstruct(struct(), categoryLabels{:});
  combo.trialFrames   = frames;
  combo.rewardDur     = [refModel.design.dspec.expt.trial.reward_duration];
  combo.baseline      = exp(full(refModel.design.X(:,refModel.design.biasCol)));
  combo.unspecialized = full(refModel.prediction);
  combo.prediction    = nan(size(refModel.prediction));
  combo.target        = refModel.target;
  combo.epoch         = refModel.design.dspec.expt.eventEpochs;
  
  %% Specialization categories and trial-specific values 
  combo.categories    = model.categories();
  combo.trialCateg    = accumfun(2, @(x) cat(1,refModel.design.dspec.expt.trial.(x)), combo.categories);
  
  %% Collect prediction of models across all trials
  for iModel = 1:numel(model)
    %% Find subset of trials predicted by this sub-model
    [~,iTrial]        = ismember(model(iModel).design.trialIndices, trialID);
    assert(all(iTrial > 0));
    
    %% Collect model prediction
    combo.prediction([frames{iTrial}])      = model(iModel).prediction;
%     [trialEpochs{iTrial}]
  end
  
  %% Canonical epoch values
  epochStart          = accumfun(2, @(x) cat(1,refModel.design.dspec.expt.trial.(x)), [cfg.behavEvents, 'duration']);
  epochDuration       = median(diff(epochStart, 1, 2), 'omitnan');
  combo.epochMSecs    = 5;
  combo.trialEpochs   = accumfun(1, @(x) butlast(linspace(x,x+1,round(epochDuration(x)/combo.epochMSecs)+1)'), 1:numel(epochDuration));

  %% Use linear interpolation to evaluate predictions at canonical epoch values
  for what = {'baseline', 'unspecialized', 'prediction'}
    prediction        = nan(numel(combo.trialEpochs), numel(frames));
    for iTrial = 1:numel(frames)
      epoch                 = combo.epoch(frames{iTrial});
      predict               = combo.(what{:})(frames{iTrial});
      sel                   = ~isnan(epoch);
      prediction(:,iTrial)  = interp1( epoch(sel), predict(sel), combo.trialEpochs, 'linear', nan );
    end
    combo.(['trial' titleCase(what{:})])    = prediction;
    
%     assert(all(isfinite(prediction(:))));
  end
  combo.trialBaseline(2:end,:)              = [];

  %% Category trial-average predictions
  combo.categoryValues= cell(size(combo.categories));
  for iCat = 1:numel(combo.categories)
    [combo.categoryValues{iCat},~,catIndex] = unique(combo.trialCateg(:,iCat));
    catLabel          = combo.categories{iCat};
    if strncmp(catLabel, 'past_', 5)
      catLabel(6)     = upper(catLabel(6));
      catLabel(5)     = [];
    end
    catLabel          = regexprep(catLabel, '_.*', '');
    isRewarded        = combo.rewardDur(:) > 0;
    
    for what = {'unspecialized', 'prediction'}
      prediction      = combo.(['trial' titleCase(what{:})]);
      catPrediction   = nan(numel(combo.trialEpochs), numel(combo.categoryValues{iCat}), 2);
      for iValue = 1:numel(combo.categoryValues{iCat})
        for iRwd = 0:1
          selTrials   = (catIndex == iValue) & (isRewarded == iRwd);
          catPrediction(:,iValue,iRwd+1)  = mean(prediction(:,selTrials) ./ combo.trialBaseline(:,selTrials), 2, 'omitnan');
        end
      end
      combo.(catLabel).(what{:})  = catPrediction;
    end
      
    %% Define modulation strength
    combo.(catLabel).modulation          ...
                      = sqrt(mean((combo.(catLabel).prediction ./ combo.(catLabel).unspecialized - 1).^2, [1 3], 'omitnan'));
%                       = mean(abs(combo.(catLabel).prediction ./ combo.(catLabel).unspecialized - 1));
  end
  
end

function iPlot = plotModelIlustration(fig, iPlot, egModel)

  %% Basis functions
  [axs,iPlot]   = fig.panel(iPlot);
  timeBin       = egModel.design.dspec.expt.binSize / 1000;   % in seconds
  time          = (0:size(egModel.design.dspec.covar(1).basis.B,1)-1) * timeBin;
  color         = hydrangea(size(egModel.design.dspec.covar(1).basis.B,2));
  
  set(axs, 'xlim', [0, egModel.design.dspec.covar(1).basis.param.duration/1000], 'ylim', [0 1.4], 'ytick', [], 'box', 'off');
  xlabel(axs, 'Time {\itt} from event {\itk} (s)');
  ylabel(axs, 'Basis {\itb_i^k}({\itt})');

  for iBasis = 1:size(egModel.design.dspec.covar(1).basis.B,2)
    line( 'parent', axs, 'xdata', time, 'ydata', egModel.design.dspec.covar(1).basis.B(:,iBasis)                                        ...
        , 'linewidth', FormatDefaults.linewidthRegular, 'color', color(iBasis,:) );
    text( egModel.design.dspec.covar(1).basis.centers(iBasis) * timeBin, 1.05, sprintf('{\\itb}_{%d}^{\\itk}',iBasis), 'color', color(iBasis,:)   ...
        , 'fontsize', FormatDefaults.legendFontSize, 'horizontalalignment', 'center', 'verticalalignment', 'bottom' );
  end

  %% Model equations
  [axs,iPlot]   = fig.panel(iPlot);
  set(axs, 'xlim', [0 2], 'ylim', [-0.3 3], 'ydir', 'reverse', 'clipping', 'off', 'xcolor', 'none', 'ycolor', 'none');

  text(0.2, -0.3, '(optional) category {\itc} of trial', 'fontsize', FormatDefaults.fontSize-1, 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'color', FormatDefaults.darkGray);
  text(0  ,  1  , 'Predict spike count {\ity}_{\itj}^{({\itc})}({\itt}) in trial {\itj}', 'fontsize', FormatDefaults.fontSize, 'horizontalalignment', 'left', 'verticalalignment', 'middle');
  text(0  ,  2  , '~ Poisson\{exp[\mu_{\itj} + \Sigma_{\itk}\Sigma_{\iti}{\itw}_{\itik}^{({\itc})}{\itb_i^k}({\itt})]\}', 'fontsize', FormatDefaults.fontSize, 'horizontalalignment', 'left', 'verticalalignment', 'middle');
%   text(0.2,  3.3, '{\itk} \in \{ fixation, C onset, delay, ... \}', 'fontsize', FormatDefaults.fontSize-1, 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'color', FormatDefaults.darkGray);
  
end

function iPlot = plotModulationStrength(fig, iPlot, model)

  cfg                 = model{1}.cfg;

  %% Average modulation strength across categorized predictions
  catLabel            = regexprep(regexprep(cfg.behavConditions, '^(past)_(.)', '$1${upper($2)}'), '_.*', '');
  modStrength         = cell(size(model));
  for iAni = 1:numel(model)
    modStrength{iAni} = nan(numel(model{iAni}.prediction), numel(catLabel));
    for iCat = 1:numel(catLabel)
      hasModulation   = arrayfun(@(x) isfield(x.(catLabel{iCat}), 'modulation'), model{iAni}.prediction);
      modStrength{iAni}(hasModulation,iCat) = arrayfun(@(x) mean(x.(catLabel{iCat}).modulation), model{iAni}.prediction(hasModulation));
    end
  end
  
  %% Distribution of modulation strengths vs. animal
  iPlot = 2;
  strengthEdges       = linspace(0, 0.8, 25);
  color               = [FormatDefaults.darkBlue; FormatDefaults.mediumPurple];
  for iCat = 1:numel(catLabel)
    [axs,iPlot]       = fig.panel(iPlot);
    
    cla(axs)
    hold(axs, 'on');
    xlabel(axs, sprintf('MI(%s)', strrep(strrep(catLabel{iCat}, 'Gap', '-gap'), 'Saccade', '-sacc')));
    if mod(iCat,2) == 1
      ylabel(axs, 'Freq. cells');
    end
    
    for iAni = 1:numel(model)
%       freq            = histcounts(modStrength{iAni}(:,iCat), strengthEdges, 'normalization', 'pdf');
%       line('parent', axs, 'xdata', strengthCenters, 'ydata', freq, 'color', color(iAni,:), 'linewidth', FormatDefaults.linewidthRegular);
      histogram(modStrength{iAni}(:,iCat), strengthEdges, 'normalization', 'prob', 'facecolor', color(iAni,:), 'edgecolor', 'none');
    end
    drawnow;
    
    %% Legend
    if iCat == 1
      axsPos          = get(axs, 'position');
      aniName         = cellfun(@(x) regexprep(x.unspecializedModel{1}.design.dspec.expt.id, '_.*', ''), model, 'UniformOutput', false);
      hLegend         = legend(axs, aniName, 'location', 'northoutside', 'numcolumns', 2);
      set(axs, 'position', axsPos);
      nudge(hLegend, [0.08 0.02]);
    end
    
    %% P-value for whether there is a significant difference across animals
    [isSigni, pValue] = kstest2(modStrength{1}(:,iCat), modStrength{2}(:,iCat));
    yRange            = get(axs, 'ylim');
    text( 0.8, yRange(2), sprintf('{\\itp} = %.1g', pValue), 'parent', axs, 'fontsize', FormatDefaults.legendFontSize   ...
        , 'horizontalalignment', 'right', 'verticalalignment', 'top', 'color', FormatDefaults.darkGray );
  end
  
end

function iPlot = plotModulationCategories(fig, iPlot, model)

  %% Get all categories by which each model is specialized
  experiment              = model{1}.hierarchicalModel{1}(1).design.dspec.expt;
  conditionCodes          = parseVariableSpecifications(experiment.desc.condition_code);

  modelSpecs              = cellfun( @(x) strjoin(formatVariableSpecifications(sort(x.categories()),conditionCodes,true), ' & ')  ...
                                   , accumfun(1, @(y) y.hierarchicalModel, model), 'UniformOutput', false);
  modelSpecs              = strrep(modelSpecs, 'past ', 'past-');
  modelSpecs              = strrep(modelSpecs, '-saccade', '-sacc');
  [modelSpecs{cellfun(@isempty, modelSpecs)}] = deal('(n.s.)');
  
  %% Replace small categories with more generic labels
  [modelType,~,typeIndex] = unique(modelSpecs);
  nOfType                 = sumDataByBin([], typeIndex, size(modelType,1))';
  
  smallFracType           = find(nOfType < 0.02*numel(modelSpecs));
  for iSmall = 1:numel(smallFracType)
    [modelSpecs{ typeIndex == smallFracType(iSmall) }]      ...
                          = deal(sprintf('misc. %d-cat.', numel(strfind(modelType{iSmall}, '&')) + 1));
  end
  
  %% Format model categories by descending order of proportion of cells (pooled across animals)
  [modelType,~,typeIndex] = unique(modelSpecs);
  nOfType                 = sumDataByBin([], typeIndex, size(modelType,1));
  [~,iOrder]              = sort(nOfType, 'descend');
  
  % Special case for small categories
  miscellaneous           = strncmp(modelType(iOrder), 'misc.', 5);
  iOrder                  = [iOrder(~miscellaneous); iOrder(miscellaneous)];  
  
  % Special case for cells w/o any significant specializations
  unspecialized           = strcmp(modelType(iOrder), '(n.s.)');
  iOrder                  = [iOrder(~unspecialized); iOrder(unspecialized)];
  
  %% Convert back to per-animal list
  typeIndex               = cellfun(@(x) typeIndex(x), span2index(cellfun(@(x) numel(x.hierarchicalModel), model)), 'UniformOutput', false);

  %% Color coding of sorted categories
  typeColor               = lines(numel(modelType));
  noPast                  = cellfun(@isempty, regexp(modelType(iOrder), '(past-|^misc.|^\(n.s.\))', 'match', 'once'));
  for iCat = find(noPast)'
    categories            = strsplit(modelType{iOrder(iCat)}, ' & ');
    pastCateg             = strcat('past-', categories);
    allCateg              = 1:numel(categories);
    pastCombo             = arrayfun(@(x) strjoin([categories(1:x-1), pastCateg(x), categories(x+1:end)], ' & '), allCateg, 'uniformoutput', false);
    pastCombo             = strrep(pastCombo, '-saccade', '-sacc');
    altColor              = pastelize(typeColor(iCat + [0,5],:), [0.3,0.7], 0);
    for iPast = 1:numel(pastCombo)
      typeColor(strcmp(modelType(iOrder), pastCombo{iPast}),:)  = altColor(iPast,:);
    end
  end
  
  miscellaneous           = strncmp(modelType(iOrder), 'misc.', 5);
  unspecialized           = strcmp(modelType(iOrder), '(n.s.)');
  typeColor(miscellaneous,:)  = repmat(linspace(0.35, 0.7, sum(miscellaneous))', 1, 3);
  typeColor(unspecialized,:)  = [1 1 1];
  
  
  %% Plot pie chart for proportions of cells with various model specializations
  for iAni = 1:numel(model)
    %% Axis formatting
    experiment            = model{iAni}.hierarchicalModel{1}(1).design.dspec.expt;
    animal                = regexprep(experiment.id, '_.*', '');
    [axs, iPlot]          = fig.panel(iPlot);
    
    %% Pie chart for proportions
    nOfType               = sumDataByBin([], typeIndex{iAni}, size(modelType,1));
    [phat, pci]           = binointerval(nOfType(iOrder), numel(typeIndex{iAni}));
    hProportion           = pie(axs, phat, unspecialized);

    colormap(axs, typeColor);
    set(axs, 'visible', 'on', 'xcolor', 'none', 'ycolor', 'none');
    hTitle                = title(axs, animal, 'fontweight', 'normal', 'horizontalalignment', 'right');
    nudge(hTitle, [-0.2 0]);
    if iAni > 1
      nudge(hTitle, [0 0.1]);
    end
    
    %% Custom labels for binomial confidence intervals on each proportion (N.B. this is not a simultaneous multinomial C.I.)
    hLabels               = findobj(hProportion, 'Type', 'text');
%     set(hLabels, 'fontsize', 13);
    
    for iType = 1:numel(phat)
      if phat(iType) < 0.03
        continue
      end
      propLabel           = sprintf( '\\fontsize{%d}%.0f_{-%.0f}^{+%.0f}%%', 12 - 2*(phat(iType) < 0.1)       ...
                                   , 100*phat(iType), 100*(phat(iType) - pci(iType,1)), 100*(pci(iType,2) - phat(iType)) );
      set(hLabels(iType), 'string', propLabel);
    end
    
    %% Legend
    if iAni > 1
      drawnow;
      axsPos              = get(axs, 'position');
      hLegend             = legend(strrep(modelType(iOrder), 'past-gap & saccade', 'past-gap & sacc'), 'location', 'northoutside', 'numcolumns', 2);
      set(axs, 'position', axsPos);
      nudge(hLegend, [0.05 0.05]);
    end
  end
  
end

function iPlot = plotModelPredictions(fig, iPlot, model)

  cfg                 = model{1}.cfg;

  %% Average modulation strength across categorized predictions
  catLabel            = regexprep(regexprep(cfg.behavConditions, '^(past)_(.)', '$1${upper($2)}'), '_.*', '');
  modStrength         = cell(size(model));
  for iAni = 1:numel(model)
    modStrength{iAni} = nan(numel(model{iAni}.prediction), numel(catLabel));
    for iCat = 1:numel(catLabel)
      hasModulation   = arrayfun(@(x) isfield(x.(catLabel{iCat}), 'modulation'), model{iAni}.prediction);
      modStrength{iAni}(hasModulation,iCat) = arrayfun(@(x) mean(x.(catLabel{iCat}).modulation), model{iAni}.prediction(hasModulation));
    end
  end
  
  %% Gap and saccade-modulated cells
  predCateg           = 'gap';
  showCateg           = cfg.behavConditions{strncmp(cfg.behavConditions, predCateg, numel(predCateg))};
  [~,modIndex]        = ismember({predCateg, 'saccade'}, catLabel);
  
  for iAni = 1:numel(model)
    %% Select cells with large modulations by all categories of interest
    modulatedCell     = find( all( modStrength{iAni}(:,modIndex) > 0, 2 )                         ...
                            & arrayfun(@(x) numel(x.trialFrames), model{iAni}.prediction)' > 80   ...
                            );
%     [~,iOrder]        = sort(sum(modStrength{iAni}(modulatedCell,modIndex), 2), 'descend');
    [~,iOrder]        = sort(modStrength{iAni}(modulatedCell,modIndex(1)), 'descend');
    modulatedCell     = modulatedCell(iOrder);

    %% Experiment info and trial categories
    experiment        = model{iAni}.unspecializedModel{1}(1).design.dspec.expt;
    animal            = regexprep(experiment.id, '_.*', '');
    catValues         = model{iAni}.prediction(modulatedCell(1)).categoryValues{strcmp(model{iAni}.prediction(modulatedCell(1)).categories, showCateg)};
    catColor          = lines(numel(catValues));

    %% Predictions for a few example cells
    for iCell = 1:3
      %% Order predictions by category values
      experiment      = model{iAni}.unspecializedModel{modulatedCell(iCell)}(1).design.dspec.expt;
      prediction      = model{iAni}.prediction(modulatedCell(iCell));
      values          = prediction.categoryValues{strcmp(prediction.categories, showCateg)};
      [~,iOrder]      = ismember(values, catValues);
    
      %% Plot prediction for correct/error trials
      for iReward = [1 0]
        [axs, iPlot]  = fig.panel(iPlot);
        
%         modulation    = prediction.(predCateg).prediction(:,iOrder,iReward+1) - prediction.(predCateg).unspecialized(:,iOrder,iReward+1);
        modulation    = prediction.(predCateg).prediction(:,iOrder,iReward+1);
        modulation    = modulation / (experiment.binSize/1000);     % convert to Hz
%         refPred       = prediction.(predCateg).unspecialized(:,iOrder,iReward+1);
%         line('parent', axs, 'xdata', prediction.trialEpochs, 'ydata', mean(refPred,2), 'linewidth', FormatDefaults.linewidthRegular, 'color', FormatDefaults.mediumGray);

        for iCat = 1:numel(iOrder)
          line('parent', axs, 'xdata', prediction.trialEpochs, 'ydata', modulation(:,iCat), 'linewidth', FormatDefaults.linewidthRegular, 'color', catColor(iCat,:));
        end
        
        %% Event onset labels
%         yRange        = [-1 1] * 1.5;
        yRange        = roundRange([0,max(modulation(:))], 0.5);
        set(axs, 'xlim', [1, numel(cfg.behavEvents)+1], 'xtick', 1:numel(cfg.behavEvents)+1, 'xgrid', 'on', 'ylim', yRange);
        if iCell == 1 && iReward == 0
          text( 0, yRange(2), 'Predicted firing rate (Hz)', 'parent', axs         ...
              , 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'fontsize', FormatDefaults.fontSize, 'rotation', 90 );
        end
        if iReward
          title(axs, sprintf('Cell %s', regexprep(experiment.id,'^[^_]*_','')), 'fontsize', FormatDefaults.legendFontSize, 'fontweight', 'normal');
          set(axs, 'xticklabel', {});
        else
          set(axs, 'xticklabel', [regexprep(regexprep(cfg.behavEvents,'_.*',''),'^c$','C onset'), 'end'], 'xticklabelrotation', 90);
        end
        
        %% Legend
        if iCell == 3 && iReward
          drawnow;
          axsPos        = get(axs, 'position');
          valueLabel    = strcat(arrayfun(@num2str, values, 'UniformOutput', false), '\circ');
          valueLabel{1} = [strrep(showCateg,'_',' '), ' = ', valueLabel{1}];
          hLegend       = legend(axs, valueLabel, 'location', 'northoutside', 'numcolumns', numel(values));
          set(axs, 'position', axsPos);
        end
      end
      
      %% Indicate error trials
      text( numel(cfg.behavEvents)+0.5, mean(yRange), 'ERROR', 'parent', axs, 'horizontalalignment', 'center', 'verticalalignment', 'middle'    ...
          , 'fontsize', FormatDefaults.legendFontSize, 'rotation', 90, 'color', FormatDefaults.mediumGray, 'fontweight', 'bold', 'background', [1 1 1] );
    end
  end
  
end
