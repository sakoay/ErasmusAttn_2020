function model = plotFiringRateDependencies(modelFile)
  
  %% Load data
  model                   = cell(size(modelFile));
  for iFile = 1:numel(modelFile)
    model{iFile}          = load(modelFile(iFile).name);
  end
  
  %% Proportions of cells with different dependencies on behavioral conditions
  plotFitDependencies(cellfun(@(x) x.hierarchicalModel, model, 'uniformoutput', false), 'GLM model dependencies');
  
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
  [pan,shape,fig]         = makePanels( numel(model), description, [description ' (proportions of cells)'], 'maxsubplotsize', 500, 'panelmargins', struct('b',5,'r',80) );
  
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
      if phat(iType) < 0.05
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
