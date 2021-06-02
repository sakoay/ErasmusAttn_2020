function model = plotFiringRatePrediction(modelFile)
  
  %% Load data
  if isstruct(modelFile)
    model             = modelFile;
  else
    model             = load(modelFile);
  end
  
  %% Proportions of cells with different dependencies on behavioral conditions
  plotFitHierarchy(model.hierarchicalModel, 'Hierarchical model', model.cfg);
  return;
  
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

function specs = formatSpecs(specs, conditions, abbreviate)
  %% Default arguments
  if ~exist('abbreviate', 'var') || isempty(abbreviate)
    abbreviate  = false;
  end

  %% Replace condition labels
  for iCond = 1:size(conditions,1)
    if abbreviate
      condStr   = regexprep(regexprep(conditions{iCond,2},' .*',''),'(rr).*','$1');
    else
      condStr   = conditions{iCond,2};
    end
    specs       = regexprep(specs, ['(condition_code=\S*)' conditions{iCond,1}], ['$1' condStr]);
  end
  
  %% Additional simplifications
  specs         = regexprep(specs, '_(direction|condition)[^=]*', '');
  specs         = regexprep(specs, '^(dir)[^=]*', '$1');
  specs         = regexprep(specs, '(_loc)ation[^=]*', '$1');
%   specs         = regexprep(specs, '^[^=]+=(.*[a-zA-Z].*)$', '$1');
  specs         = regexprep(specs, '^condition_code=', '');
  specs         = regexprep(specs, '_+', ' ');
end

function fig = plotPredictionByCondition(model, description, cfg, showDetails)

  %% Default arguments
  if ~exist('showDetails', 'var') || isempty(showDetails)
    showDetails           = true;
  end
  
  %% Use the basic unspecialized model as a comparison
  refModel                = unique(model.firstAncestor);
  experiment              = refModel.design.dspec.expt;
%   yRange                  = [0, max(refModel.target)];
  yRange                  = [0, quantile(refModel.target,0.999)];
  rateTicks               = 1:ceil((yRange(2) - 1)/4):yRange(2);
  
  %% Specify categories of trials for which to plot data together
  trialConditions         = accumfun(2, @(x) cat(1,experiment.trial.(x)), cfg.behavConditions);
  [trialCategory, ~, categIndex]    = unique(trialConditions, 'rows');
  
  %% Labels for trial categories
  categLabel              = arrayfun(@(x) parseVariableSpecifications(experiment.desc.(cfg.behavConditions{x}),trialCategory(:,x)), 1:numel(cfg.behavConditions), 'UniformOutput', false);
  categLabel              = [categLabel{:}];
  
  %% Configure plots
  [pan,shape,fig]         = makePanels( size(trialCategory,1) + 1, experiment.id, [strrep(experiment.id,'_',' ') ' : ' description]     ...
                                      , 'aspectratio', 1, 'plotmargins', struct('b',20), 'panelmargins', struct('b',15,'t',15), 'maxsubplotsize', 200 );

  dataColor               = [1 1 1] * 0.5;
  if showDetails
    refColor              = [1 1 1] * 0.1;
    modelColor            = [0 0 0];
    basisColor            = [ 0.1373    0.6471    0.3333    ... fix_on
                            ; 0         0.4470    0.7410    ... cue_start
                            ; 0.3010    0.7450    0.9330    ... cue_stop
                            ; 0.8500    0.3250    0.0980    ... c_start
                            ; 0.6350    0.0780    0.1840    ... c_stop
                            ; 0.9290    0.6940    0.1250    ... delay_start
                            ; 0.4940    0.1840    0.5560    ... choice_target
                            ; 0.6706    0.3529    0.6667    ... saccade_onset
                            ; 0.4660    0.6740    0.1880    ... reward_start
                            ];
  else
    refColor              = [0 123 255] / 255;
    modelColor            = [0.8 0 0];
  end
  
  %% Loop through categories of trials and plot each one
  for iCat = 1:size(trialCategory,1)
    %% Select time bins corresponding to all trials within this category
    trialIndices          = find(categIndex == iCat);
    trial                 = experiment.trial(trialIndices);
    time                  = 0:max([trial.duration]);
    
    trialBins             = experiment.trialBins(trialIndices);
    selData               = [trialBins{:}];
    nTrials               = numel(trialBins);
    nTimeBins             = max(cellfun(@numel, trialBins));
    timeBin               = toBinCenters( (0:nTimeBins) * experiment.binSize );
    dataIdx               = accumfun(2, @(x) nTimeBins*(x-1) + 1 + trialBins{x} - trialBins{x}(1), 1:nTrials);

    %% Axes formatting
    axs                   = selectPanel(pan, iCat, shape);
    hold(axs, 'on');
    axis(axs, 'tight');
    colormap(axs, flipud(gray));
    title(axs, sprintf('%d %s : %s', numel(trial), categLabel{iCat,1}, strjoin(categLabel(iCat,2:end),' ')), 'fontweight', 'normal');
    ylabel(axs, sprintf('# spikes / %.3g%s', experiment.binSize, experiment.unitOfTime));
    set(axs, 'xcolor', [1 1 1]*0.5, 'ycolor', [1 1 1]*0.5, 'clim', [0 1], 'clipping', 'on');

    %% Use grid to indicate time scale
    xGrid               = 0:1000:timeBin(end);           % 1 second intervals
    line( 'parent', axs, 'xdata', colvec(repmat(xGrid,3,1)), 'ydata', colvec(repmat([yRange(:); nan],1,numel(xGrid))), 'linestyle', '-.', 'color', [1 1 1]*0.7 );

    %% Plot mean firing rate across trials of this category
    firingRate            = nan(nTimeBins, numel(trialBins));
    firingRate(dataIdx)   = refModel.target(selData);
    
    bandplot(axs, timeBin, firingRate, [], dataColor);
%     line( 'parent', axs, 'xdata', colvec(timeBin(:) + 0.1*experiment.binSize*randn(1,nTrials))            ...
%         , 'ydata', colvec(firingRate + 0.2*randn(1,nTrials)), 'color', [1 1 1]*0.7, 'marker', '.', 'markersize', 10, 'linestyle', 'none');
%     line('parent', axs, 'xdata', timeBin, 'ydata', mean(firingRate,2,'omitnan'), 'color', [1 1 1]*0.3, 'linewidth', 1);
    
    %% Plot reference vs. model predictions for this category
    refPredict            = nan(nTimeBins, numel(trialBins));
    prediction            = refModel.jointPrediction();
    refPredict(dataIdx)   = prediction(selData);
    modelPredict          = nan(nTimeBins, numel(trialBins));
    prediction            = model.jointPrediction();
    modelPredict(dataIdx) = prediction(selData);

    line('parent', axs, 'xdata', timeBin, 'ydata', mean(modelPredict,2,'omitnan'), 'color', modelColor, 'linewidth', 2);
    line('parent', axs, 'xdata', timeBin, 'ydata', mean(refPredict,2,'omitnan')  , 'color', refColor  , 'linewidth', 2, 'linestyle', '--');
    
    %% Indicate stimulus periods with a line
    yOffset               = yRange(end) * 0.05;
    iOffset               = 0;
    ticks                 = [];
    tickLabels            = {};
    for what = {'cue', 'c'}
      onTime              = zeros(numel(time), numel(trial));
      for iTrial = 1:numel(trial)
        if isnan(trial(iTrial).([what{:} '_start']))
          continue
        end
        onTime( 1 + (trial(iTrial).([what{:} '_start']):trial(iTrial).([what{:} '_stop'])), iTrial )  = 1;
      end
      onTime              = mean(onTime, 2);

      iOffset             = iOffset + 1;
      ticks(end+1)        = -iOffset * yOffset;
      tickLabels{end+1}   = what{:};
      image('parent', axs, 'xdata', time([1 end]), 'ydata', -(iOffset + 0.15*[-1 1])*yOffset, 'cdata', onTime', 'CDataMapping', 'scaled');
    end
    
    set( axs, 'ytick', [flip(ticks), rateTicks], 'yticklabel', [flip(tickLabels), arrayfun(@num2str,rateTicks,'uniformoutput',false)]   ...
       , 'ylim', [-0.1 1]*yRange(end), 'XAxisLocation', 'origin' );
%        , 'ylim', [-0.3 1]*yRange(end), 'XAxisLocation', 'origin' );

    %% Indicate various behavioral events
    behavEvent          = {'fix_on', 'delay_start', 'choice_target', 'reward_start'};
%     eventColor          = lines(numel(behavEvent));
    ticks               = [];
    tickLabels          = {};
    for iEvent = 1:numel(behavEvent)
      ticks(end+1)      = median([trial.(behavEvent{iEvent})], 'omitnan');
%       tickLabels{end+1} = coloredText(regexprep(behavEvent{iEvent}, '_.+', ''), eventColor(iEvent,:));
      tickLabels{end+1} = regexprep(behavEvent{iEvent}, '_.+', '');
    end

    if any(isnan(ticks))
      xlabel(axs, 'Time (ms)');
    else
      [ticks,iOrder]    = sort(ticks);
      set(axs, 'xtick', ticks, 'xticklabel', tickLabels(iOrder), 'xticklabelrotation', 90, 'xlim', [0 timeBin(end)]);
    end
    
    %% Identify the single model that is used to predict this trial category
    if ~showDetails
      continue;
    end
    
    catModel              = model( arrayfun(@(x) all(ismember(trialIndices,x.design.trialIndices)), model) );
    assert(numel(catModel) == 1);

    %% Plot regressors aligned to the median event time associated with the regressor type
    eventName             = fieldnames(catModel.regressors);
    eventName(~cellfun(@(x) isstruct(catModel.regressors.(x)), eventName))  = [];
    
%     sumResponses          = zeros(numel(timeBin),1);
    for iReg = 1:numel(eventName)
      eventTime           = median([trial.(eventName{iReg})], 'omitnan');     % N.B. some events may not be defined for a particular trial, in which case the regressor does not apply to the prediction for that trial
      regressor           = catModel.regressors.(eventName{iReg});
      regBin              = experiment.binfun(eventTime) + (0:numel(regressor.time)-1);
      sel                 = regBin <= numel(timeBin);
      regBin(~sel)        = [];
      
      %% Plot the contribution to the Poisson distribution mean that is solely due to this regressor
      % N.B. This assumes that we performed a GLM fit with Poisson distributed data and a
      % logarithmic link function
%       line('parent', axs, 'xdata', eventTime + regressor.time, 'ydata', exp(regressor.response), 'color', basisColor(iReg,:), 'linewidth', 1);
      line('parent', axs, 'xdata', timeBin(regBin), 'ydata', exp(catModel.regressors.bias + regressor.response(sel)), 'color', basisColor(iReg,:), 'linewidth', 1);
%       sumResponses(regBin)= sumResponses(regBin) + regressor.response(sel);
    end
%     line('parent', axs, 'xdata', timeBin, 'ydata', exp(catModel.regressors.bias + sumResponses), 'color', [240 200 46]/255, 'linewidth', 1, 'linestyle', '-.');
  end
  
  %% Legend
  axs                   = selectPanel(pan, size(trialCategory,1) + 1, shape);
  set(axs, 'xcolor', 'none', 'ycolor', 'none');
  hPlot                 = gobjects(0);
  hPlot(end+1)          = line('parent', axs, 'xdata', [], 'ydata', [], 'linewidth', 1, 'color', dataColor);
  hPlot(end+1)          = line('parent', axs, 'xdata', [], 'ydata', [], 'linewidth', 2, 'color', refColor, 'linestyle', '--');
  hPlot(end+1)          = line('parent', axs, 'xdata', [], 'ydata', [], 'linewidth', 2, 'color', modelColor);
  for iReg = 1:numel(eventName)
    hPlot(end+1)        = line('parent', axs, 'xdata', [], 'ydata', [], 'color', basisColor(iReg,:), 'linewidth', 1);
  end
  columnlegend(2, [{'Data \pm std.dev.'; 'Unspecialized model'; 'Specialized model'}; strrep(eventName,'_',' ')], 'object', hPlot, 'location', 'east');
  
end

function fig = plotRegressorsByCondition(model, description, cfg, showDetails)

  %% Default arguments
  if ~exist('showDetails', 'var') || isempty(showDetails)
    showDetails           = false;
  end
  
  %% Use the basic unspecialized model as a comparison
  refModel                = unique(model.firstAncestor);
  experiment              = refModel.design.dspec.expt;
  conditionCodes          = parseVariableSpecifications(experiment.desc.condition_code);

  eventName               = fieldnames(refModel.regressors);
  eventName(~cellfun(@(x) isstruct(refModel.regressors.(x)), eventName))  = [];

  yRange                  = cellfun(@(x) max(exp(refModel.regressors.bias + refModel.regressors.(x).response)), eventName);

  %% Specify categories of trials for which to plot regressors together
  modelSpecs              = model.specializations(true);
  categories              = fieldnames(modelSpecs);
  [catValues,~,catIndex]  = cellfun(@(x) unique({modelSpecs.(x)}), categories, 'UniformOutput', false);
  nSubCat                 = cellfun(@numel, catValues(2:end));
  if isempty(nSubCat)
    nSubCat               = 1;
    subCatIndex           = ones(size(model));
  else
    subCatIndex           = sub2ind([nSubCat; 1], catIndex{2:end});
  end
  
  yRange(:,2:prod(nSubCat)) = 0;
  
  %% Configure plots
  lineStyle               = {'-.', '-'};
  catColor                = lines(numel(catValues{1}));
%   if strcmp(categories{1}, 'condition_code')
%     [~,iCond]             = ismember(catValues{1}, conditionCodes(:,1));
%     catValues{1}          = conditionCodes(iCond,2)';
%   end
  catLabel                = parseVariableSpecifications(experiment.desc.(categories{1}), catValues{1});
  catLabel                = [strrep(categories{1},'_',' '), ' = ', strjoin(coloredText(catLabel, catColor), ' | ')];
  
  [pan,shape,fig]         = makePanels( [prod(nSubCat),numel(eventName)], experiment.id, [strrep(experiment.id,'_',' ') ' : ' catLabel]       ...
                                      , 'aspectratio', 1, 'plotmargins', struct('b',25), 'panelmargins', struct('b',15,'t',20), 'maxsubplotsize', 180 );
  
  %% Compare each regressor across trial categories
  for iReg = 1:numel(eventName)
    for iCat = 1:prod(nSubCat)
      iModel              = find(subCatIndex == iCat,1,'first');
      if isempty(iModel)
        continue
      end
      
      %% Category label
      specs               = model(iModel).specializations(true);
      catLabel            = cellfun(@(x) specs.(x), categories(2:end), 'UniformOutput', false);
      catLabel            = cellfun(@(x,y) parseVariableSpecifications(experiment.desc.(x),y,true), categories(2:end), catLabel, 'UniformOutput', false);
      sel                 = cellfun(@iscell, catLabel);
      catLabel(sel)       = cellfun(@(x) ['[' strjoin(x,';') ']'], catLabel(sel), 'UniformOutput', false);
      catLabel            = strcat(categories(2:end), '=', catLabel);
      catLabel            = strjoin(formatSpecs(catLabel, conditionCodes), ', ');

      %% Axes formatting
      axs                 = selectPanel(pan, [iCat,iReg], shape);
      hold(axs, 'on');
      axis(axs, 'tight');
      title(axs, catLabel, 'fontweight', 'normal');
%       xlabel(axs, sprintf('t from %s (%s)', strrep(eventName{iReg},'_',' '), experiment.unitOfTime));
      xlabel(axs, sprintf('t from %s', strrep(eventName{iReg},'_',' ')));
      ylabel(axs, sprintf('# spikes / %.3g%s', experiment.binSize, experiment.unitOfTime));

      %% Plot reference unspecialized model regressors
      regressor           = refModel.regressors.(eventName{iReg});
      line( 'parent', axs, 'xdata', regressor.time, 'ydata', exp(model(iModel).regressors.bias + regressor.response)  ...
          , 'color', [1 1 1]*0.7, 'linewidth', 2, 'linestyle', '--' );
      
      %% Plot the contribution to the Poisson distribution mean that is solely due to this regressor
      % N.B. This assumes that we performed a GLM fit with Poisson distributed data and a
      % logarithmic link function
      for iModel = 1:numel(model)
        regressor         = model(iModel).regressors.(eventName{iReg});
        selModel          = subCatIndex(iModel) == iCat;
        if ~showDetails && ~selModel
          continue
        end
        
        response          = exp(model(iModel).regressors.bias + regressor.response);
        yRange(iReg,iCat) = max([yRange(iReg,iCat); response]);
        line( 'parent', axs, 'xdata', regressor.time, 'ydata', response   ...
            , 'color', catColor(catIndex{1}(iModel),:), 'linewidth', 1+selModel, 'linestyle', lineStyle{selModel+1} );
      end
    end
  end
  
  %% Synchronize plot ranges
  isOutlier               = yRange > 5*median(yRange(:));
  yRange                  = [0, max(yRange(~isOutlier))];

  for iReg = 1:numel(eventName)
    for iCat = 1:prod(nSubCat)
      if isOutlier(iReg,iCat)
        continue
      end
      set(selectPanel(pan, [iCat,iReg], shape), 'ylim', yRange);
    end
  end
  
end

function fig = plotFitHierarchy(model, description, cfg)

  %% Get hierarchy descriptions for all models in the given set
  model(cellfun(@isempty, model)) = [];
  modelSpecs              = cell(numel(model), 0);
  for iModel = 1:numel(model)
    specs                 = model{iModel}.hierarchy();
    modelSpecs(iModel,1:numel(specs))  = specs;
  end
  
  %% UGLY : condense category labels
  experiment              = model{1}(1).design.dspec.expt;
  conditionCodes          = parseVariableSpecifications(experiment.desc.condition_code);
  modelSpecs              = formatSpecs(modelSpecs, conditionCodes, true);
  
  %% Configure plots
  animal                  = regexprep(experiment.id, '_.*', '');
  [pan,shape,fig]         = makePanels( 2, animal, [animal ' : ' description], 'maxsubplotsize', 1000, 'aspectratio', 2 );
  set(fig, 'Tag', '<persist>');
  
  %% Plot detailed vs. condensed hierarchy tree
  warning('off', 'MATLAB:table:RowsAddedExistingVars');
  for iDetail = 0:1
    axs                   = selectPanel(pan, iDetail + 1, shape);
    if iDetail
      specialization      = modelSpecs;
      nudge(axs, [-0.1 0 0.15 0]);
    else
      specialization      = regexprep(modelSpecs, '=.*', '');
      nudge(axs, [-0.05 0 0 0]);
    end
    modelFits             = model;
    
    %% Define directed graph with edge weights being the number of cells with the same sequence of category specializations
    parent                = repmat({animal}, size(specialization,1), 1);
    modelGraph            = digraph(0, table(parent(1),{model},{{}}, 'VariableNames', {'Name','Model','Plot'}));
    for iDepth = 1:size(specialization,2)
      %% Identify trivial nodes at this depth (cells with no further specializations)
      notSpecialized      = cellfun(@isempty, specialization(:,iDepth));
      [specialization{notSpecialized,iDepth}] = deal('');

      %% Construct node names as the entire hierarchy up to this depth of specialization
      nodes               = strcat(parent, {'->'}, specialization(:,iDepth));
      [nodes,cellIdx,nodeIdx] = unique(nodes);
      nodeModels          = arrayfun(@(x) modelFits(nodeIdx==x), 1:numel(nodes), 'UniformOutput', false);
      nCells              = cellfun(@numel, nodeModels);

      %% Add edges with weights proportional to the number of cells with the corresponding type of specialization
      modelGraph          = modelGraph.addedge(parent(cellIdx), nodes, nCells);
      [~,iNode]           = ismember(nodes, modelGraph.Nodes.Name);
      modelGraph.Nodes(iNode,2)         = nodeModels(:);

      %% Update hierarchy information up to this depth
      parent              = nodes(nodeIdx);
      parent(notSpecialized)            = [];
      specialization(notSpecialized,:)  = [];
      modelFits(notSpecialized)         = [];
    end
    
    %% Prune trivial branches for simplicity
    % Find all leaf nodes that correspond to no further specialization
    leafIdx               = find(modelGraph.outdegree() == 0);
    leafIdx(~endsWith(modelGraph.Nodes{leafIdx,'Name'}, '->'))  = [];

    % Prune trivial leaf nodes where its parent node only has it as a child
    doPrune               = [];
    for iLeaf = 1:numel(leafIdx)
      [~,parentIdx]       = modelGraph.inedges(leafIdx(iLeaf));
      if modelGraph.outdegree(parentIdx) < 2
        doPrune(end+1)    = leafIdx(iLeaf);
      end
    end
    modelGraph            = modelGraph.rmnode(doPrune);

    %% Construct and display graph object
    set(axs, 'xcolor', 'none', 'ycolor', 'none');
    axis(axs, 'tight');
    
    hGraph                = plot(modelGraph, 'Parent', axs, 'Layout', 'layered', 'LineWidth', 0.5*modelGraph.Edges.Weight);
    % Condense node labels
    nodeLabel             = get(hGraph, 'NodeLabel');
    labeledge(hGraph, 1:size(modelGraph.Edges,1), modelGraph.Edges.Weight);
    set( hGraph, 'PickableParts', 'none', 'NodeLabel', regexprep(nodeLabel, '^.*->([^>]*)$', '$1'), 'NodeFontSize', 12, 'MarkerSize', 7   ...
       , 'EdgeColor', [1 1 1]*0.6, 'EdgeLabelColor', [1 1 1]*0.3, 'NodeLabelColor', [0.8 0 0], 'UserData', modelGraph );
    
    %% Interactive plotting of cells with the selected hierarchy
    hSelect               = line('Parent', axs, 'XData', [], 'YData', [], 'Marker', 'o', 'MarkerSize', 11, 'LineWidth', 3, 'MarkerEdgeColor', [240 200 46]/255, 'PickableParts', 'none');
    set(axs, 'Box', 'off', 'XTick', [], 'YTick', [], 'ButtonDownFcn', {@showCellsInHierarchy,hGraph,hSelect,description,cfg}, 'Interruptible', 'off', 'BusyAction', 'cancel');
  end
  warning('on', 'MATLAB:table:RowsAddedExistingVars');
  
end

function fig = compareModelPredictions(model, description, cfg)
  %% Configure plots
  [pan,shape,fig]         = makePanels( size(trialCategory,1), experiment.id, [strrep(experiment.id,'_',' ') ' : ' description]     ...
                                      , 'aspectratio', 1, 'plotmargins', struct('b',20), 'panelmargins', struct('b',15,'t',15), 'maxsubplotsize', 200 );
  refColor                = [0 123 255] / 255;
  modelColor              = [0.8 0 0];
end

function showCellsInHierarchy(hObject, event, hGraph, hSelect, varargin)
  
  %% Locate the closest graph node to where the user has clicked
  clickPos              = get(hObject, 'CurrentPoint');
  clickPos              = clickPos(1,1:2);
  nodePos               = [get(hGraph,'XData'); get(hGraph,'YData')]';
  dist2                 = sum(bsxfun(@minus, nodePos, clickPos).^2,2);
  [minDist2,iClosest]   = min(dist2);
  
  if minDist2 > 0.5^2
    set(hSelect, 'XData', [], 'YData', []);
    return;
  end
  
  %% Retrieve models corresponding to the selected node
  nodePos               = nodePos(iClosest,:);
  set(hSelect, 'XData', nodePos(1), 'YData', nodePos(2));
  
  modelGraph            = get(hGraph, 'UserData');
  models                = modelGraph.Nodes{iClosest, 'Model'}{:};
  
  %% Don't remake the plot if it already exists
  hFig                  = modelGraph.Nodes{iClosest, 'Plot'}{:};
  if ~isempty(hFig) && all(ishghandle(hFig(:)))
    for iModel = 1:numel(models)
      figure(hFig(iModel));
    end
  else
    if ~isempty(hFig)
      delete(hFig);
    end
    set(get(hObject,'Parent'), 'Pointer', 'watch');
    drawnow;
    
    hFig                = gobjects(numel(models),2);
    for iModel = 1:numel(models)
%       hFig(iModel,1)    = plotPredictionByCondition(models{iModel}, varargin{:});
      hFig(iModel,2)    = plotRegressorsByCondition(models{iModel}, varargin{:});
      drawnow;
    end
    
    modelGraph.Nodes{iClosest, 'Plot'}  = {hFig};
    set(hGraph, 'UserData', modelGraph);
    set(get(hObject,'Parent'), 'Pointer', 'arrow');
  end

  assignin('base', 'selModel', models);
  
end
