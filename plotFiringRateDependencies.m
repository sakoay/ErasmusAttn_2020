function model = plotFiringRateDependencies(modelFile)
  
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

  %% Proportions of cells with different dependencies on behavioral conditions
  hierarchicalModel         = cellfun(@(x) x.hierarchicalModel, model, 'uniformoutput', false);
  plotFitDependencies(hierarchicalModel, 'GLM model dependencies');
  
  %% 
  plotCategoryPredictions(hierarchicalModel, cellfun(@(x) x.prediction, model, 'uniformoutput', false), model{1}.cfg);
  
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

function fig = plotCategoryPredictions(model, combo, cfg)


  %% Category trial-average predictions
  catLabel            = regexprep(regexprep(cfg.behavConditions, '^(past)_(.)', '$1${upper($2)}'), '_.*', '');
  modStrength         = cell(size(combo));
  for iAni = 1:numel(combo)
%     modStrength{iAni} = -0.1 * ones(numel(combo{iAni}), numel(catLabel));
    modStrength{iAni} = nan(numel(combo{iAni}), numel(catLabel));
    for iCat = 1:numel(catLabel)
      hasModulation   = arrayfun(@(x) isfield(x.(catLabel{iCat}), 'modulation'), combo{iAni});
      modStrength{iAni}(hasModulation,iCat) = arrayfun(@(x) mean(x.(catLabel{iCat}).modulation), combo{iAni}(hasModulation));
    end
  end
  
  %%
  for iCat = 1:numel(catLabel)
    figure; hold on; 
    for iAni = 1:numel(combo)
      histogram(min(modStrength{iAni}(:,iCat),1,'includenan'), linspace(0,1,50), 'normalization', 'pdf');
    end
    [isSigni, pValue] = kstest2(modStrength{1}(:,iCat), modStrength{2}(:,iCat));
    title(sprintf('%s : %.3g', catLabel{iCat}, pValue));
  end
  
  %%
  for iAni = 1:numel(combo)
    saccadeCell       = find(modStrength{iAni}(:,1) > 0 & modStrength{iAni}(:,end) > 0);
    animal            = regexprep(model{iAni}{1}(1).design.dspec.expt.id, '_.*', '');
    [pan,shape]       = makePanels(numel(saccadeCell), animal, [animal ': gap & saccade modulated']);
    for iCell = 1:numel(saccadeCell)
      axs             = selectPanel(pan,iCell,shape);
      plot(axs, combo{iAni}(saccadeCell(iCell)).trialEpochs, combo{iAni}(saccadeCell(iCell)).gap.prediction(:,:,end) - combo{iAni}(saccadeCell(iCell)).gap.unspecialized(:,:,end));
%       hold(axs, 'on');
%       plot(axs, combo{iAni}(saccadeCell(iCell)).trialEpochs,combo{iAni}(saccadeCell(iCell)).gap.prediction(:,:,end));
%       plot(axs, combo{iAni}(saccadeCell(iCell)).trialEpochs, combo{iAni}(saccadeCell(iCell)).gap.unspecialized(:,:,end), 'k');
      set(axs, 'xlim', [1, numel(cfg.behavEvents)+1], 'xtick', 1:numel(cfg.behavEvents)+1, 'xticklabel', [regexprep(cfg.behavEvents,'_.*',''), 'end'], 'xticklabelrotation', 90, 'ylim', [-1 1]);
      ylabel(axs, 'specialized - unspecialized');
      title(axs, mean(combo{iAni}(saccadeCell(iCell)).(catLabel{iCat}).modulation));
    end
  end
  
end

function fig = plotFitDependencies(model, description)

  for iAni = 1:numel(model)
    model{iAni}(cellfun(@isempty, model{iAni})) = [];
  end

  %% Get all categories by which each model is specialized
  experiment              = model{1}{1}(1).design.dspec.expt;
  conditionCodes          = parseVariableSpecifications(experiment.desc.condition_code);

  modelSpecs              = cellfun(@(x) strjoin(formatVariableSpecifications(sort(x.categories()),conditionCodes,true), ' & '), cat(1,model{:}), 'UniformOutput', false);
  modelSpecs              = strrep(modelSpecs, 'past ', 'past-');
  modelSpecs              = strrep(modelSpecs, '-saccade', '-sacc');
  [modelSpecs{cellfun(@isempty, modelSpecs)}] = deal('(n.s.)');
  
  %% Replace small categories with more generic labels
  [modelType,~,typeIndex] = unique(modelSpecs);
  nOfType                 = sumDataByBin([], typeIndex, size(modelType,1))';
  
  smallFracType           = find(nOfType < 0.015*numel(modelSpecs));
  for iSmall = 1:numel(smallFracType)
    [modelSpecs{ typeIndex == smallFracType(iSmall) }]      ...
                          = deal(sprintf('mixed(%d)', numel(strfind(modelType{iSmall}, '&')) + 1));
  end
  
  %% Format model categories by descending order of proportion of cells (pooled across animals)
  [modelType,~,typeIndex] = unique(modelSpecs);
  nOfType                 = sumDataByBin([], typeIndex, size(modelType,1));
  [~,iOrder]              = sort(nOfType, 'descend');
  
  % Special case for small categories
  miscellaneous           = strncmp(modelType(iOrder), 'mixed', 5);
  iOrder                  = [iOrder(~miscellaneous); iOrder(miscellaneous)];  
  
  % Special case for cells w/o any significant specializations
  unspecialized           = strcmp(modelType(iOrder), '(n.s.)');
  iOrder                  = [iOrder(~unspecialized); iOrder(unspecialized)];
  
  %% Color coding of sorted categories
  typeColor               = lines(numel(modelType));
  typeColor(8:end,:)      = pastelize(typeColor(8:end,:), [0.3,0.7], 0.3);
  miscellaneous           = strncmp(modelType(iOrder), 'mixed', 5);
  unspecialized           = strcmp(modelType(iOrder), '(n.s.)');
  typeColor(miscellaneous,:)  = repmat(linspace(0.3, 0.7, sum(miscellaneous))', 1, 3);
  typeColor(unspecialized,:)  = [1 1 1];
  
  %% Convert back to per-animal list
  typeIndex               = cellfun(@(x) typeIndex(x), span2index(cellfun(@numel, model)), 'UniformOutput', false);

  
  %% Configure plots
  [pan,shape,fig]         = makePanels( numel(model), description, [description ' (proportions of cells)'], 'maxsubplotsize', 500, 'panelmargins', struct('b',5,'r',100) );
  
  %% Plot pie chart for proportions of cells with various model specializations
  for iAni = 1:numel(model)
    %% Axis formatting
    experiment            = model{iAni}{1}(1).design.dspec.expt;
    animal                = regexprep(experiment.id, '_.*', '');
    axs                   = selectPanel(pan, iAni, shape);
    title(axs, animal);
    colormap(axs, typeColor);
    
    %% Pie chart for proportions
    nOfType               = sumDataByBin([], typeIndex{iAni}, size(modelType,1));
    [phat, pci]           = binointerval(nOfType(iOrder), numel(typeIndex{iAni}));
    hProportion           = pie(axs, phat, unspecialized);
    
    %% Custom labels for binomial confidence intervals on each proportion (N.B. this is not a simultaneous multinomial C.I.)
    hLabels               = findobj(hProportion, 'Type', 'text');
%     set(hLabels, 'fontsize', 13);
    
    for iType = 1:numel(phat)
      if phat(iType) < 0.03
        continue
      end
      propLabel           = sprintf('\\fontsize{13}%.0f_{-%.0f}^{+%.0f}%%', 100*phat(iType), 100*(phat(iType) - pci(iType,1)), 100*(pci(iType,2) - phat(iType)));
      set(hLabels(iType), 'string', propLabel);
    end
    
    %% Legend
    if iAni > 1
      hLegend             = legend(modelType(iOrder), 'location', 'eastoutside');
      nudge(axs, [-0.15 0 0.2 0]);
      nudge(hLegend, [0.04 0]);
    else
      nudge(axs, [-0.05 0]);
    end
  end
    
end
