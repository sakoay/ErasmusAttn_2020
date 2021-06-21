function fitBehaviorData(dataFile, postfix)

  %% Default arguments
  if ~exist('postfix', 'var') || isempty(postfix)
    postfix           = '';
  end

  %% Analysis configuration
  cfg                 = struct();

  % Mapping of gap direction onto -1/+1
  cfg.gapSignMap      = [  270     -1     ... left
                        ;   90     +1     ... right
                        ;  180     -1     ... down
                        ;    0     +1     ... up
                        ];
  %                         -1     +1         sign mapping
  cfg.ClocationMap    = [    5      1     ... left/right -> horizontal
                        ;    3      7     ... down/up    -> vertical
                        ];
  
  % Cross validation and permutation tests
  cfg.nBootstraps     = 100;
  cfg.nCVFolds        = 5;
  cfg.nLambdas        = 20;
  
  % Trial selection criteria
  cfg.selectConditions= 1:2;            % keep only trials with these condition_code
  cfg.numHistoryTerms = 3;              % number of past-trial conditions to include in model
  
  % Configuration for lassoglm()
  cfg.target          = 'condition_code';
  cfg.fitOptions      = {'binomial', 'Link', 'logit', 'MaxIter', 1e5, 'Alpha', 0.5};
  
  %% Load data
  if isstruct(dataFile)
    data              = dataFile;
  else
    data              = load(dataFile);
    data              = data.behavior;
  end
  
  %% Define output file
  [path,name,ext]     = parsePath(dataFile);
  name                = sprintf('behavModel_%s_%dpast', regexprep(name, '^[^_]+_', ''), cfg.numHistoryTerms);
  outputFile          = fullfile(path, [name postfix ext]);
  
  
  %% Configure models per cell
  experiment          = cell(numel(data), 1);
  Xvariables          = [ {'gapLR'}, arrayfun(@(x) sprintf('past%d_gapLR',x), 1:cfg.numHistoryTerms, 'UniformOutput', false)  ...
                        , {'gapUD'}, arrayfun(@(x) sprintf('past%d_gapUD',x), 1:cfg.numHistoryTerms, 'UniformOutput', false)  ...
                        , arrayfun(@(x) sprintf('past%d_saccadeLR',x), 1:cfg.numHistoryTerms, 'UniformOutput', false)         ...
                        , arrayfun(@(x) sprintf('past%d_saccadeUD',x), 1:cfg.numHistoryTerms, 'UniformOutput', false)         ...
                        ];
  targetY             = cell(numel(data), 2);
  behaviorX           = cell(numel(data), 2);
  nLocations          = size(cfg.ClocationMap,1);
  parfor iCell = 1:numel(data)
    %% Select data for a single cell
    experiment{iCell} = [data(iCell).monkey '_' data(iCell).cell_id];
    
    %% Apply trial selection to cell data, including the requirement that trials have the desired number of historical info
    selTrials         = ismember(data(iCell).condition_code.n, cfg.selectConditions);
    if ~any(selTrials)
      continue;
    end

    %% Separately collect data for horizontal vs. vertical gap direction possibilities
    for iLoc = 1:nLocations
      %% Define target as the saccade direction
      sel                     = selTrials & ismember(data(iCell).direction_C.n, cfg.ClocationMap(iLoc,:));
      behavX                  = nan(sum(sel), numel(Xvariables));
      targetY{iCell,iLoc}     = data(iCell).saccade_direction.n(sel);
      iX                      = 0;
      
      %% Collect gap directions up to the desired number of history terms
      for iDir = 1:nLocations
        iX                    = iX + 1;
        behavX(:,iX)          = data(iCell).gap_direction.n(sel);
        for iPast = 1:cfg.numHistoryTerms
          iX                  = iX + 1;
          behavX(:,iX)        = data(iCell).gap_direction.(sprintf('n_%d',iPast))(sel);
        end
        
        %% Flag with NaNs gap directions that don't match the iDir selection (horizontal vs. vertical)
        for jX = iX - cfg.numHistoryTerms:iX
          [~,gapIndex]        = ismember(behavX(:,jX), cfg.gapSignMap(:,1));
          behavX( ~ismember(gapIndex, 2*iDir + (-1:0)), jX )  = nan;
        end
      end
      
      %% Collect saccade directions up to the desired number of history terms
      for iDir = 1:nLocations
        for iPast = 1:cfg.numHistoryTerms
          iX                  = iX + 1;
          behavX(:,iX)        = data(iCell).saccade_direction.(sprintf('n_%d',iPast))(sel);
        end
        
        %% Flag with NaNs saccade directions that don't match the iDir selection (horizontal vs. vertical)
        for jX = iX - cfg.numHistoryTerms + 1:iX
          [~,gapIndex]        = ismember(behavX(:,jX), cfg.gapSignMap(:,1));
          behavX( ~ismember(gapIndex, 2*iDir + (-1:0)), jX )  = nan;
        end
      end
      
      %% Flag with NaNs past trial info where the gap direction was of a different modality than the current trial)
      %{
      for iPast = 1:cfg.numHistoryTerms
        isDiffMode            = ~ismember(data(iCell).direction_C.(sprintf('n_%d',iPast))(sel), cfg.ClocationMap(iLoc,:));
        behavX(isDiffMode,[0,cfg.numHistoryTerms]+1+iPast)  = nan;
      end
      %}
      
      %% Convert orientation info in degrees to a -1/+1 map
      targetY{iCell,iLoc}     = arrayfun(@(x) zeroIfEmpty(cfg.gapSignMap(x==cfg.gapSignMap(:,1),2)), targetY{iCell,iLoc});
      behaviorX{iCell,iLoc}   = arrayfun(@(x) zeroIfEmpty(cfg.gapSignMap(x==cfg.gapSignMap(:,1),2)), behavX);
    end
  end
  
  %% Pool data across cells for each monkey
  for name = unique({data.monkey})
    selData           = strcmp({data.monkey}, name{:});
    experiment{end+1} = name{:};
    targetY(end+1,:)  = catcell(1, targetY(selData,:), 1);
    behaviorX(end+1,:)= catcell(1, behaviorX(selData,:), 1);
  end
  
  %% Logistic regression
  behaviorModel       = cell(size(targetY));
  parfor iData = 1:numel(behaviorModel)
    %% Require a minimum number of trials for fitting the model
    if numel(targetY{iData}) < cfg.nCVFolds * numel(Xvariables)
      continue
    end
    
    %% Use cross-validation to determine regularization strength
    X                 = behaviorX{iData};
    Y                 = targetY{iData} > 0;
    [modelW,fitInfo]  = lassoglm(X, Y, cfg.fitOptions{:}, 'NumLambda', cfg.nLambdas, 'CV', cfg.nCVFolds);
    lambdaRatio       = fitInfo.LambdaMinDeviance / fitInfo.Lambda(end);
%     lassoPlot(modelW,fitInfo,'plottype','CV');

    %% Use bootstrap experiments to obtain uncertainties on fitted weights
    [~,bootstrapExp]  = permutationExperiments(numel(Y), 0, cfg.nBootstraps, cfg.nCVFolds);
    bootstrapExp(:,1) = [];               % omit original experiment
    if lambdaRatio < 1
      bootstrapW      = nan(size(modelW,1), cfg.nBootstraps);
      for iBS = 1:cfg.nBootstraps
        bsW           = lassoglm(X(bootstrapExp(:,iBS),:), Y(bootstrapExp(:,iBS)), cfg.fitOptions{:}, 'NumLambda', 2, 'LambdaRatio', lambdaRatio);
        bootstrapW(:,iBS) = bsW(:,1);
      end
    else
      bootstrapW      = zeros(size(modelW,1), cfg.nBootstraps);
    end
    
    %% Collect fit information
    fitInfo.Xvariables= Xvariables;
    fitInfo.designX   = X;
    fitInfo.targetY   = Y;
    fitInfo.modelW    = modelW;
    fitInfo.bootstrapW= bootstrapW;
    behaviorModel{iData}  = fitInfo;
    
%     longfigure; imagesc(bootstrapW); tempscale
  end
  
  %% Save output
  save(outputFile, 'experiment', 'behaviorModel', 'cfg');
  fprintf(' -->  %s\n', outputFile);
  
end

function x = zeroIfEmpty(x)
  if isempty(x)
    x = 0;
  end
end
