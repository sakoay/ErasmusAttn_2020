function model = plotEncodingModes(modelFile)
  
  %% Load data
  if isstruct(modelFile)
    model             = modelFile;
  elseif iscell(modelFile)
    model             = cellfun(@load, modelFile);
  else
    model             = load(modelFile);
  end
  
  %% Analysis configuration
  cfg                 = model.cfg;
  cfg.maxDepth        = 2;
  cfg.minNCells       = 5;
  
  %% Proportions of cells with different dependencies on behavioral conditions
  plotEventRegressors(cat(1,model.hierarchicalModel), cfg);
%   plotSpecializationStrengths(cat(1,model.hierarchicalModel), cfg);
%   plotRegressorModes(cat(1,model.hierarchicalModel), cfg);
  
end

function fig = plotEventRegressors(model, cfg)

  model(cellfun(@isempty, model)) = [];
  
  %% Identify all specializations up to depth 2 for which there are a minimum number of cells with this specialization
  specCategories          = cellfun(@(x) x.categories(), model, 'UniformOutput', false);
  specializations         = {};
  specCells               = {};
  for iDepth = 0:cfg.maxDepth
    %% Include only cells with at least depth of specialization
    subSpecs              = cellfun(@(x) x(1:min(iDepth,end)), specCategories, 'UniformOutput', false);
    subSpecs(cellfun(@numel, subSpecs) < iDepth)  = [];
    subSpecs              = cellfun(@sort, subSpecs, 'UniformOutput', false);     % order of categories don't matter
    
    %% Identify cells with relevant hierarchies
    [~,egIndex,specIndex] = unique(cellfun(@(x) strjoin(x,' & '), subSpecs, 'UniformOutput', false));
    selSpecs              = find( arrayfun(@(x) sum(specIndex==x), 1:numel(egIndex)) >= cfg.minNCells );
    specializations       = [specializations; subSpecs(egIndex(selSpecs))];
    specCells             = [specCells; arrayfun(@(x) find(specIndex==x), selSpecs(:), 'UniformOutput', false)];
  end
  
  %% 
  for iSpec = 1:numel(specializations)
    %% Pool experimental data across models to determine all valid categories of trials
    specs                 = specializations{iSpec};
    if isempty(specs)
      specs               = {''};
      conditions          = nan;
    else
      trialConditions     = accumfun(1, @(x) accumfun(2, @(y) cat(1,x(1).design.dspec.expt.trial.(y)), specs), model);
      conditions          = unique(trialConditions, 'rows');
    end
    
    %% Get regressors specialized for each of the category values in this set
    [response,time]       = accumfun(1, @(x) x.responseByCategory(specs,conditions), model(specCells{iSpec}));
    behavEvents           = fieldnames(response);
    
    %% Define which conditions to overlay vs. make separate plots for
    [compCond,~,condIdx]  = unique(conditions(:,1));
    [otherCond,~,otherIdx]= unique(conditions(:,2:end), 'rows');
    condLabel             = formatVariableSpecifications(specs);

    %% Plot separately for other conditions
    for iOther = 1:size(otherCond,1)
      %% Configure plots
      otherLabel          = arrayfun(@(x) sprintf('%s = %.3g',condLabel{x},otherCond(iOther,x-1)), 2:numel(condLabel), 'UniformOutput', false);
      [pan,shape,fig]     = makePanels( [numel(compCond), numel(behavEvents)], strjoin(specs,'_'), strjoin(otherLabel,' & '), 'maxsubplotsize', 150 );
      
      %% Plot regressors for each behavioral event
      for iEvent = 1:numel(behavEvents)
        for iCond = 1:numel(compCond)
          %%
          axs             = selectPanel(pan, [iCond,iEvent]);
          if ~isempty(specs{1})
            title(axs, sprintf('%s = %.3g',condLabel{1},compCond(iCond)), 'FontWeight', 'normal');
          end
          xlabel(axs, ['t from ' strrep(behavEvents{iEvent},'_',' ')]);
          ylabel(axs, 'Response');
        
          %%
          eventResponse   = accumfun(2, @(x) x.(behavEvents{iEvent})(:,condIdx==iCond & otherIdx==iOther), response);
          eventTime       = accumfun(2, @(x) time.(behavEvents{iEvent}), time);
          eventResponse(end+1,:)  = nan;
          eventTime(end+1,:)      = nan;
          patch('parent', axs, 'xdata', eventTime, 'ydata', eventResponse, 'facecolor', 'none', 'edgecolor', [0.8 0 0], 'linewidth', 1, 'edgealpha', 0.5);

          line('parent', axs, 'xdata', eventTime([1 end-1],1), 'ydata', [0 0], 'linestyle', '-.', 'color', [1 1 1]*0.7, 'linewidth', 2);
        end
      end
      drawnow;
    end
  end
  
end

function fig = plotSpecializationStrengths(model, cfg)

  model(cellfun(@isempty, model)) = [];
  
  %% Identify all specializations up to depth 2 for which there are a minimum number of cells with this specialization
  specCategories          = cellfun(@(x) x.categories(), model, 'UniformOutput', false);
  specializations         = {};
  specCells               = {};
  for iDepth = 0:cfg.maxDepth
    %% Include only cells with at least depth of specialization
    subSpecs              = cellfun(@(x) x(1:min(iDepth,end)), specCategories, 'UniformOutput', false);
    subSpecs(cellfun(@numel, subSpecs) < iDepth)  = [];
    subSpecs              = cellfun(@sort, subSpecs, 'UniformOutput', false);     % order of categories don't matter
    
    %% Identify cells with relevant hierarchies
    [~,egIndex,specIndex] = unique(cellfun(@(x) strjoin(x,' & '), subSpecs, 'UniformOutput', false));
    selSpecs              = find( arrayfun(@(x) sum(specIndex==x), 1:numel(egIndex)) >= cfg.minNCells );
    specializations       = [specializations; subSpecs(egIndex(selSpecs))];
    specCells             = [specCells; arrayfun(@(x) find(specIndex==x), selSpecs(:), 'UniformOutput', false)];
  end
  
  %% Population modes for major forms of specialization
  for iSpec = 1:numel(specializations)
    %% Pool experimental data across models to determine all valid categories of trials
    specs                 = specializations{iSpec};
    if isempty(specs)
      specs               = {''};
      conditions          = nan;
    else
      trialConditions     = accumfun(1, @(x) accumfun(2, @(y) cat(1,x(1).design.dspec.expt.trial.(y)), specs), model);
      conditions          = unique(trialConditions, 'rows');
    end
    
    %% Get regressors specialized for each of the category values in this set
    [response,time]       = accumfun(1, @(x) x.responseByCategory(specs,conditions), model(specCells{iSpec}));
    behavEvents           = fieldnames(response);
%     cellResponse          = accumfun(2, @(x) remould(permute(cat(3,response.(x)),[3 1 2]),2), behavEvents);
    cellResponse          = accumfun(2, @(x) sqsum(cat(3,response.(x)),1)', behavEvents);

    %% 
    cellDist              = squareform(pdist(cellResponse, 'euclidean'));
    clusterTree           = linkage(cellDist, 'single');
    cellOrder             = optimalleaforder(clusterTree, cellDist);
%     figure;   imagesc(cellDist(cellOrder,cellOrder));   title(strrep(strjoin(specs,' & '),'_',' '));
    figure;   imagesc(cellResponse(cellOrder,:));  tempscale;   title(strrep(strjoin(specs,' & '),'_',' '));

    
    %% Configure plots
%     animal                = regexprep(experiment.id, '_.*', '');
%     [pan,shape,fig]       = makePanels( 2, animal, [animal ' : ' description], 'maxsubplotsize', 1000, 'aspectratio', 2 );

  end
  tilefigs
  
end

function fig = plotRegressorModes(model, cfg)

  model(cellfun(@isempty, model)) = [];
  
  %% Identify all specializations up to depth 2 for which there are a minimum number of cells with this specialization
  specCategories          = cellfun(@(x) x.categories(), model, 'UniformOutput', false);
  specializations         = {};
  specCells               = {};
  for iDepth = 0:cfg.maxDepth
    %% Include only cells with at least depth of specialization
    subSpecs              = cellfun(@(x) x(1:min(iDepth,end)), specCategories, 'UniformOutput', false);
    subSpecs(cellfun(@numel, subSpecs) < iDepth)  = [];
    subSpecs              = cellfun(@sort, subSpecs, 'UniformOutput', false);     % order of categories don't matter
    
    %% Identify cells with relevant hierarchies
    [~,egIndex,specIndex] = unique(cellfun(@(x) strjoin(x,' & '), subSpecs, 'UniformOutput', false));
    selSpecs              = find( arrayfun(@(x) sum(specIndex==x), 1:numel(egIndex)) >= cfg.minNCells );
    specializations       = [specializations; subSpecs(egIndex(selSpecs))];
    specCells             = [specCells; arrayfun(@(x) find(specIndex==x), selSpecs(:), 'UniformOutput', false)];
  end
  
  %% Population modes for major forms of specialization
  for iSpec = 1:numel(specializations)
    %% Pool experimental data across models to determine all valid categories of trials
    specs                 = specializations{iSpec};
    if isempty(specs)
      specs               = {''};
      conditions          = nan;
    else
      trialConditions     = accumfun(1, @(x) accumfun(2, @(y) cat(1,x(1).design.dspec.expt.trial.(y)), specs), model);
      conditions          = unique(trialConditions, 'rows');
    end
    
    %% Get regressors specialized for each of the category values in this set
    [response,time]       = accumfun(1, @(x) x.responseByCategory(specs,conditions), model(specCells{iSpec}));
    behavEvents           = fieldnames(response);
%     cellResponse          = accumfun(2, @(x) remould(permute(cat(3,response.(x)),[3 1 2]),2), behavEvents);
    cellResponse          = accumfun(2, @(x) sqsum(cat(3,response.(x)),1)', behavEvents);

    %% 
    cellDist              = squareform(pdist(cellResponse, 'euclidean'));
    clusterTree           = linkage(cellDist, 'single');
    cellOrder             = optimalleaforder(clusterTree, cellDist);
%     figure;   imagesc(cellDist(cellOrder,cellOrder));   title(strrep(strjoin(specs,' & '),'_',' '));
    figure;   imagesc(cellResponse(cellOrder,:));  tempscale;   title(strrep(strjoin(specs,' & '),'_',' '));

    
    %% Configure plots
%     animal                = regexprep(experiment.id, '_.*', '');
%     [pan,shape,fig]       = makePanels( 2, animal, [animal ' : ' description], 'maxsubplotsize', 1000, 'aspectratio', 2 );

  end
  tilefigs
  
end
