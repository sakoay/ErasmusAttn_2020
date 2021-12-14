% Regression model for a subset of trials specified by some category selection criteria
classdef CategorizedModel < handle
  
  %_________________________________________________________________________________________________
  properties (SetAccess = protected)
    ID                = []            % unique identifier for each model, not changed upon copying
    
    progenitor        = []            % the parent model that includes data across all values of this model's category
    categoryLabel     = ''            % the trial category label that is used to select data used in this model fit
    categoryValues    = []            % values of categoryLabel that are included in the selected data in the fit for this model

    target            = []            % the time series data predicted by this model
    design            = []            % design matrix specifying covariates used to predict target 
    prediction        = []            % cross-validated prediction for target 
    fitInfo           = []            % fit information for model, with cross-validated goodness of fit (fitInfo.Deviance)
    regressors        = struct()      % fitted weights, time axis and waveforms for regressors

    cvTrainSel        = []            % selection vector for subsets of target used to generate training sets for cross-validation
    extra             = struct()      % optional user specified information
  end
  
  %_________________________________________________________________________________________________
  methods (Static)
    
    %----- Apply additional trial selections to design matrix
    function [design, selData] = selectTrials(design, selTrials)
      design.trialIndices(~selTrials)       = [];

      %% Select timepoints corresponding to the selected trials
      selData                               = ismember(design.trialX, design.trialIndices);
      design.X(~selData,:)                  = [];
      design.trialX(~selData,:)             = [];

      %% Remove constant columns of the new design matrix
      if ~isempty(design.X)
        isConst                             = full(all(design.X(1,:) == design.X(2:end,:), 1));
        design.constCols(~design.constCols) = isConst;
        design.X(:,isConst)                 = [];
      end
    end
    
    %----- Structure conversion to load an object of this class from disk
    function obj = loadobj(frozen)
      %% HACK : I have no idea why this can sometimes already be the target type
      if isa(frozen, 'CategorizedModel')
        obj = frozen;
        return;
      end

      %% Start from default constructor 
      obj(size(frozen,1), size(frozen,2))  = CategorizedModel();

      %% Inherit all frozen fields 
      for field = fieldnames(frozen)'
        for iObj = 1:numel(frozen)
          obj(iObj).(field{:})  = frozen(iObj).(field{:});
        end
      end
      
      %% HACK : ensure that category labels are strings
      for iObj = 1:numel(frozen)
        if isempty(obj(iObj).categoryLabel)
          obj(iObj).categoryLabel = '';
        end
      end
      
      %% Restore raw data from progenitor models
      for iObj = 1:numel(obj)
        obj(iObj).ID      = CategorizedModel.nextID();
        if ~obj(iObj).hasProgenitor()
          continue;
        end
        
        [obj(iObj).design, selData]                                                                                                 ...
                          = CategorizedModel.selectTrials( obj(iObj).progenitor.design                                              ...
                                                         , ismember(obj(iObj).progenitor.design.trialIndices, frozen(iObj).design)  ...
                                                         );
        obj(iObj).target  = obj(iObj).progenitor.target(selData);
      end
    end
    
    %----- Get the next available ID for unique model identification
    function id = nextID(reset)
      persistent tally;
      if nargin > 0
        tally   = reset;
      elseif isempty(tally)
        tally   = 1;
      else
        tally   = tally + 1;
      end
      id        = tally;
    end

    
    %----- Log likelihood for Poisson distributed data
    function lnL = poissonLogLikelihood(target, prediction, varargin)
      lnL = -prediction + target .* log(prediction) - gammaln(target + 1);
      lnL = sum(lnL, 1);
    end
    
    %----- Deviance for Poisson distributed data: twice the difference between the log-likelihood of the fitted model and the saturated model
    function dev = poissonDeviance(target, prediction, varargin)
      dev = 2*( target .* log((target + (target==0)) ./ prediction) - (target - prediction) );
      dev = sum(dev, 1);
    end
    
    %----- Bayesian Information Criteria for LASSO
    function bic = poissonBIC(target, prediction, activeSet)
      neg2LnL       = -CategorizedModel.poissonLogLikelihood(target, prediction);
      if islogical(activeSet)
        dof         = sum(activeSet, 1);
      else
        dof         = activeSet;
      end
      bic           = neg2LnL + dof * (log(size(target,1)) / size(target,1));
    end
    
    %----- Aikaike Information Criteria for LASSO
    function aic = poissonAIC(target, prediction, activeSet)
      neg2LnL       = -CategorizedModel.poissonLogLikelihood(target, prediction);
      if islogical(activeSet)
        dof         = sum(activeSet, 1);
      else
        dof         = activeSet;
      end
      aic           = neg2LnL + dof * (2 / size(target,1));
    end
  
    
    %----- Optimize regularized log-likelihood for a series of regularization strengths 
    function [modelW, fitInfo] = fitRegularizationPath(y, X, lambdas, options)
      modelW                        = nan(size(X,2), numel(lambdas));
      W0                            = zeros(size(X,2), 1);
      fitInfo                       = struct('cost', nan(size(lambdas)), 'lambda', lambdas, 'iterations', nan(size(lambdas)));
      for iLambda = numel(lambdas):-1:1
        options.lambda              = lambdas(iLambda);
        if nargout > 1
          [modelW(:,iLambda),optim] = mexFistaFlat(y, X, W0, options);
          fitInfo.cost(iLambda)     = optim(1,:);
          fitInfo.iterations(iLambda) = optim(4,:);
        else
          modelW(:,iLambda)         = mexFistaFlat(y, X, W0, options);
        end
        
        %% Warm start for next lambda in the list
        W0                          = modelW(:,iLambda);
      end
    end

    %----- Cross-validated (unpenalized) negative log-likelihood for a series of regularization strengths 
    function [negLnL, modelW, prediction, lambdaMax] = xvalRegularizationPath(y, X, lambdaRatio, options, cvTrainSel)
      %% Loop over cross-validation folds
      modelW                        = nan(size(X,2), numel(lambdaRatio), size(cvTrainSel,2));
      prediction                    = nan(numel(y), numel(lambdaRatio));
      lambdaMax                     = nan(1, size(cvTrainSel,2));
      for iCV = 1:size(cvTrainSel,2)
        trainSel                    = cvTrainSel(:,iCV);
        
        %% Align regularization strengths across folds by adjusting to the per-fold guesstimated maximum
        selData                     = trainSel & y > 0;
        lambdaMax(iCV)              = max( X(selData,1:end-options.intercept)' * ( log(y(selData)) - X(selData,end-options.intercept+1:end) ) );
        lambdas                     = lambdaRatio * lambdaMax(iCV);
        
        %% Train model on a subset of data
        modelW(:,:,iCV)             = CategorizedModel.fitRegularizationPath(y(trainSel), X(trainSel,:), lambdas, options);

        %% Evaluate model on left-out test data
        testSel                     = ~trainSel;
        prediction(testSel,:)       = exp(X(testSel,:) * modelW(:,:,iCV));
      end
      
      %% Negative log-likelihood
      negLnL                        = -CategorizedModel.poissonLogLikelihood(y, prediction);
    end
    
  end
    
  %_________________________________________________________________________________________________
  methods
    
    %----- Constructor with no model components
    function obj = CategorizedModel(target, design)
      if nargin < 1
        return;
      end
      
      obj.ID          = CategorizedModel.nextID();
      obj.target      = target;
      if isfield(design, 'constCols')
        obj.design    = design;
      else
        obj.design    = buildGLM.removeConstantCols(design);
      end
      if size(obj.target,1) ~= size(design.X,1)
        error('CategorizedModel:constructor', 'target data must have as many rows as the design matrix, was %d vs. %d instead.', size(obj.target,1), size(design.X,1));
      end
    end
    
    %----- Structure version to store an object of this class to disk
    function frozen = saveobj(this)
      %% Use class metadata to determine what properties to save
      metadata    = metaclass(this);
      
      %% Store all mutable and non-transient data
      frozen      = repmat(struct, size(this));
      for iProp = 1:numel(metadata.PropertyList)
        property  = metadata.PropertyList(iProp);
        if property.Transient || property.Constant
          continue;
        end
        
        for iObj = 1:numel(this)
          frozen(iObj).(property.Name)  = this(iObj).(property.Name);
        end
      end
      
      %% Omit storage of raw data when it can be recreated from a progenitor model
      for iObj = 1:numel(this)
        if this(iObj).hasProgenitor()
          frozen(iObj).target = [];
          frozen(iObj).design = this(iObj).design.trialIndices;
        end
      end
    end
    
    %----- Locate unique models by ID
    function [C,ia,ic] = unique(obj)
      [~, ia, ic] = unique([obj.ID]);
      C           = obj(ia);
    end
    
    %----- Display information and inheritance tree
    function print(obj, prefix)
      %% Default arguments
      if exist('prefix', 'var') && ~isempty(prefix)
        label         = prefix;
      else
        prefix        = '';
      end

      %% Print each model in this set
      for iObj = 1:numel(obj)
        if isempty(prefix)
          label       = sprintf('%3d -- ', iObj);
        end
        
        fprintf('%s[%5d]  %*s%4d trials : ', label, obj(iObj).ID, 20-numel(label), '', obj(iObj).numberOfTrials());
        if isempty(obj(iObj).categoryLabel)
          fprintf('%33s', '');
        else
          fprintf('%15s = %-15s', obj(iObj).categoryLabel, mat2str(obj(iObj).categoryValues'));
        end
        if ~isempty(obj(iObj).fitInfo)
          fprintf(' ... %3d parameters, lnL = %.10g', obj(iObj).effectiveDegreesOfFreedom(), obj(iObj).goodnessOfFit());
        end
        fprintf('\n');
        if obj(iObj).hasProgenitor()
          obj(iObj).progenitor.print([repmat(' ',1,conditional(isempty(prefix),5,numel(prefix))), '\--']);
        end
      end
    end

    
    %----- Whether this model is a specialization of a more general model that includes trials with all values of categoryLabel
    function yes = hasProgenitor(obj)
      yes           = false(size(obj));
      for iObj = 1:numel(obj)
        if ~isempty(obj(iObj).progenitor)
          yes(iObj) = true;
        end
      end
    end
    
    %----- Get the list of unique IDs for all parent models of the given set, or nan for models that have none
    function id = parentID(obj)
      id            = nan(size(obj));
      for iObj = 1:numel(obj)
        if ~isempty(obj(iObj).progenitor)
          id(iObj)  = obj(iObj).progenitor.ID;
        end
      end
    end
    
    %----- Follow the inheritance chain upwards to the first nonzero model
    function ancestor = firstAncestor(this)
      %% Loop over all models 
      ancestor            = this;
      for iObj = 1:numel(this)
        %% Recursively crawl up the inheritance chain
        parent            = this(iObj).progenitor;
        while ~isempty(parent)
          ancestor(iObj)  = parent;
          parent          = parent.progenitor;
        end
      end
    end
    
    %----- Get the number of trials in the data specified for this model
    function count = numberOfTrials(obj)
      count           = nan(size(obj));
      for iObj = 1:numel(obj)
        if ~isempty(obj(iObj).design)
          count(iObj) = numel(obj(iObj).design.trialIndices);
        end
      end
    end
    
    %----- Get the effective number of degrees of freedom i.e. nonzero L1-regularized coefficients
    function dof = effectiveDegreesOfFreedom(obj)
      dof           = nan(size(obj));
      for iObj = 1:numel(obj)
        if isempty(obj(iObj).fitInfo)
        elseif isfield(obj(iObj).fitInfo, 'modelW')
          dof(iObj) = sum(obj(iObj).fitInfo.modelW ~= 0, 1);
        else
          dof(iObj) = sum(obj(iObj).fitInfo.activeSet, 1);
        end
      end
    end
    
    %----- Get the values of all specializations that define this model 
    function specs = specializations(obj, asString)
      %% Default arguments
      if ~exist('asString', 'var') || isempty(asString)
        asString      = false;
      end

      %% Check that all objects have the same category specialization
      if any(~strcmp(obj(1).categoryLabel, {obj(2:end).categoryLabel}))
        error('CategorizedModel:specializations', 'Inconsistent category specializations for objects in this set.');
      end
      
      %% Inherit specializations of parent
      if obj(1).hasProgenitor
        specs         = obj(1).progenitor.specializations(asString);
        for iObj = 2:numel(obj)
          specs(iObj) = obj(iObj).progenitor.specializations(asString);
        end
        specs         = reshape(specs, size(obj));
      else
        specs         = repmat(struct(), size(obj));
      end
      
      %% Register category values
      if ~isempty(obj(1).categoryLabel)
        for iObj = 1:numel(obj)
          if asString
            specs(iObj).(obj(iObj).categoryLabel) = mat2str(obj(iObj).categoryValues);
          else
            specs(iObj).(obj(iObj).categoryLabel) = obj(iObj).categoryValues;
          end
        end
      end
    end
    
    %----- Get the set of all trial categories that this model is specialized for
    function labels = categories(obj)
      %% Check that all objects have the same category specialization
      if any(~strcmp(obj(1).categoryLabel, {obj(2:end).categoryLabel}))
        error('CategorizedModel:categories', 'Inconsistent category specializations for objects in this set.');
      end
      
      if isempty(obj(1).categoryLabel)
        labels  = {};
      else
        labels  = {obj(1).categoryLabel};
      end
      
      %% Concatenate categories for parent models
      if obj(1).hasProgenitor()
        parent  = [obj.progenitor];
        labels  = [parent.categories(), labels];
      end
    end
    
    %----- Describe this set of models as the hierarchy of category specializations
    function specs = hierarchy(obj)
      if obj(1).hasProgenitor()
        %% Check that all objects have the same category specialization
        if any(~strcmp(obj(1).categoryLabel, {obj(2:end).categoryLabel}))
          error('CategorizedModel:hierarchy', 'Inconsistent category specializations for objects in this set.');
        end
        
        %% Check that the category values across objects form a partition of possible values (no overlaps)
        [catSets,objIdx]  = unique(arrayfun(@(x) strrep(mat2str(rowvec(x.categoryValues)),' ',','), obj, 'UniformOutput', false));
        allValues         = cat(1, obj(objIdx).categoryValues);
        if numel(unique(allValues)) ~= numel(allValues)
          error('CategorizedModel:hierarchy', 'Specialization values do not form a partition of %s.', obj(1).categoryLabel);
        end
        
        %% Concatenate hierarchy descriptor for the parent model
        parent            = unique([obj.progenitor]);
        parentSpecs       = parent.hierarchy();
        specs             = [parentSpecs, {[obj(1).categoryLabel, '=', strjoin(catSets,'|')]}];
      elseif numel(obj) > 1
        error('CategorizedModel:hierarchy', 'Invalid number %d of models with no parent.', numel(obj));
      elseif isempty(obj.categoryLabel)
        specs                 = '';
      else
        specs                 = [obj.categoryLabel '=', mat2str(obj.categoryValues)];
      end
    end

    
    %----- Create more specialized models from all possible partitions of the data by values of a given category
    function subModels = specialize(obj, categoryLabel, allPartitions, minNumTrials)
      %% Default arguments
      if ~exist('allPartitions', 'var') || isempty(allPartitions)
        allPartitions               = false;
      end
      if ~exist('minNumTrials', 'var') || isempty(minNumTrials)
        minNumTrials                = 5;
      end
      
      %% All models are required to be based on the same data
      experiment                    = obj(1).design.dspec.expt;
      if any(~arrayfun(@(x) isequal(x.design.dspec.expt.id,experiment.id), obj(2:end)))
        error('CategorizedModel:specialize', 'Inconsistent experimental data encountered for the given set of models.');
      end
      
      %% Identify all partitions of values of dataCategory in the dataset for this model
      dataCategory                  = cat(1, experiment.trial.(categoryLabel));
      category                      = unique(dataCategory(accumfun(2, @(x) x.design.trialIndices, obj)));
      if allPartitions
        categorySets                = partitions(category);
      else
        categorySets                = {num2cell(category)};
      end
      categorySets( cellfun(@numel, categorySets) < 2 ) = [];     % disallow trivial partitions

      %% Create models per partition
      subModels                     = repmat({repmat(CategorizedModel,0)}, size(categorySets));
      for iSet = 1:numel(categorySets)
        for iObj = 1:numel(obj)      
          for iVal = 1:numel(categorySets{iSet})
            %% Apply further trial selection to design matrix
            [design, selData]       = CategorizedModel.selectTrials(obj(iObj).design, ismember(dataCategory(obj(iObj).design.trialIndices), categorySets{iSet}{iVal}));
            if isempty(design.X)    % no data in this category
              continue;
            end

            %% Construct a more specialized model if possible
            model                   = CategorizedModel(obj(iObj).target(selData), design);
            model.progenitor        = obj(iObj);
            model.categoryLabel     = categoryLabel;
            model.categoryValues    = categorySets{iSet}{iVal};
            if model.numberOfTrials() < minNumTrials || model.numberOfTrials() == obj(iObj).numberOfTrials()
              %% If there is insufficient data for model fitting, or no specialization, inherit the parent model for this subset of data
              model.prediction      = obj(iObj).prediction(ismember(obj(iObj).design.trialX, design.trialIndices));
              model.fitInfo         = obj(iObj).fitInfo;
              model.regressors      = obj(iObj).regressors;
            end
            subModels{iSet}(end+1)  = model;
          end
        end
      end
    end
    
    %----- Fit regression model
    function obj = regress(obj, nCVFolds, fitOptions, lazy, nLambdas, minLambda)
      %% Default arguments
      if ~exist('lazy', 'var') || isempty(lazy)
        lazy                          = false;
      end
      if ~exist('nLambdas', 'var') || isempty(nLambdas)
        nLambdas                      = 10;
      end
      if ~exist('minLambda', 'var') || isempty(minLambda)
        minLambda                     = 1e-4;
      end

      %% Fit model for each configuration of model specifications and data to explain
      cvTrainSel                      = cell(size(obj));
      fitInfo                         = cell(size(obj));
      modelW                          = cell(size(obj));
      cvPrediction                    = cell(size(obj));
      parfor iObj = 1:numel(obj)
        if lazy && ~isempty(obj(iObj).prediction)
          continue;
        end
        
        %% Generate lambda ratio sequence on a log scale 
        lambdaRatio                   = exp(linspace(log(minLambda), log(1), nLambdas));

        %% GLM configuration
        options                       = fitOptions;
        options.lambda                = lambdaRatio;
        options.groups                = int32(accumfun(2, @(x) x*ones(1,obj(iObj).design.dspec.covar(x).edim), 1:numel(obj(iObj).design.dspec.covar)));
%         options.groups(end+1:size(obj(iObj).design.X,2))  = 0;
        if ~isempty(obj(iObj).design.biasCol)
          options.groups(obj(iObj).design.biasCol)        = 0;
        end
        options.groups(obj(iObj).design.constCols)        = [];
      
        %% Define cross-validation selection and shuffled data
        trialLength                   = SplitVec(obj(iObj).design.trialX, 'equal', 'length');
        [~,~,cvTrainSel{iObj}]        = permutationExperiments(trialLength, 0, 0, nCVFolds, 2);
        xvalPartitions                = cvpartition(numel(trialLength), 'KFold', nCVFolds);
        designX                       = obj(iObj).design.X;
        targetY                       = obj(iObj).target;
        
%         fitglm(designX, targetY, 'linear', 'Distribution', 'poisson', 'Link', 'log', 'Intercept', true );
%         fitInfo                   = linearSupportVectorMR( obj(iObj).target, obj(iObj).design.X, obj(iObj).cvTrainSel, obj(iObj).shuffleExp, bootstrapExp, [], [], [], true );

        %% Use cross-validation to determine regularization strength
        fitInfo{iObj}.cv1NegLnL       = CategorizedModel.xvalRegularizationPath(targetY, designX, lambdaRatio, options, cvTrainSel{iObj}(:,:,1));
        [~,bestIndex]                 = min(fitInfo{iObj}.cv1NegLnL);

        %% Use all of data to obtain central value for model weights
        % Guesstimate range of regularization strengths
        selData                       = obj(iObj).target > 0;
        lambdaMax                     = max( obj(iObj).design.X(selData,1:end-options.intercept)'                                           ...
                                           * ( log(obj(iObj).target(selData)) - obj(iObj).design.X(selData,end-options.intercept+1:end) )   ...
                                           );
        [modelW{iObj},fitInfo{iObj}]  = CategorizedModel.fitRegularizationPath(targetY, designX, lambdaRatio(bestIndex) * lambdaMax, options);
        fitInfo{iObj}.bestIndex       = bestIndex;
        
%         longfigure;   hold on;    plot(targetY);  plot(exp(designX * modelW{iObj}));    title(fitInfo{iObj}.cost)
%         figure; plot(modelW{iObj}'); figure; plot(fitInfo{iObj}.cv1NegLnL);

        %% Use an independent set of cross-validation folds to evaluate goodness of fit
        [fitInfo{iObj}.cv2NegLnL, fitInfo{iObj}.cvW, cvPrediction{iObj}]          ...
                                      = CategorizedModel.xvalRegularizationPath(targetY, designX, lambdaRatio(bestIndex), options, cvTrainSel{iObj}(:,:,2));
        fitInfo{iObj}.cvW             = squish(fitInfo{iObj}.cvW, 2);
      end
      
      %% Store fitted parameters
      for iObj = 1:numel(obj)
        if lazy && isempty(fitInfo{iObj})
          continue;
        end
        
        obj(iObj).fitInfo             = fitInfo{iObj};
        obj(iObj).fitInfo.modelW      = modelW{iObj};
        obj(iObj).prediction          = cvPrediction{iObj};
        obj(iObj).regressors          = buildGLM.combineWeights(obj(iObj).design, modelW{iObj});
%         obj(iObj).regressors.bias     = modelW{iObj}(end - options.intercept + 1:end);
      end
    end
    
    %----- Combined prediction of the given set of category-specialized models
    function [predict, allTarget] = jointPrediction(obj)
      %% Special case of a single model
      if numel(obj) == 1
        predict           = obj.prediction;
        if nargout > 1
          allTarget       = obj.target;
        end
        return;
      end
      
      %% All models are required to be based on the same data
      experiment          = obj(1).design.dspec.expt;
      if any(~arrayfun(@(x) isequal(x.design.dspec.expt.id,experiment.id), obj(2:end)))
        error('CategorizedModel:jointPrediction', 'Inconsistent experimental data encountered for the given set of models.');
      end

      %% Not necessary : The union of all models in this set are required to explain all data
      %{
      nTrials             = numel(obj(1).design.dspec.expt.trial);
      trials              = accumfun(1, @(x) x.trialIndices(:), [obj.design]);
      if numel(trials) ~= nTrials || ~isequal(sort(trials), (1:nTrials)')
        error('CategorizedModel:jointPrediction', 'Union of models in the given set must explain all of the data (%d trials).', nTrials);
      end
      %}
      
      %% Collect predictions from all models
      predict             = nan(experiment.numTimeBins, 1);
      if nargout > 1
        allTarget         = nan(experiment.numTimeBins, 1);
      end
      
      for iObj = 1:numel(obj)
        bins              = [obj(iObj).design.dspec.expt.trialBins{obj(iObj).design.trialIndices}];
        predict(bins)     = obj(iObj).prediction;
        if nargout > 1
          allTarget(bins) = obj(iObj).target;
        end
      end
    end
    
    %----- Cross-validated goodness of fit: gofFcn(target,prediction,activeSet) can be poissonLogLikelihood() %etc.
    function gof = goodnessOfFit(obj, gofFcn)
      %% Default arguments
      if ~exist('gofFcn', 'var') || isempty(gofFcn)
        gofFcn            = @CategorizedModel.poissonLogLikelihood;
      end
      
      %% Preallocate output 
      nExperiments        = max(arrayfun(@(x) size(x.prediction,2), obj));
      if nExperiments > 1
        gof               = nan(nExperiments, numel(obj));
      else
        gof               = nan(size(obj));
      end
      
      %% Loop over models
      for iObj = 1:numel(obj)
        if isempty(obj(iObj).prediction)
          continue;
        end
        if nExperiments > 1
          gof(:,iObj)     = gofFcn(obj(iObj).target, obj(iObj).prediction, obj(iObj).effectiveDegreesOfFreedom());
        else
          gof(iObj)       = gofFcn(obj(iObj).target, obj(iObj).prediction, obj(iObj).effectiveDegreesOfFreedom());
        end
      end
    end
    
    %----- Relative likelihood of this set of category-specialized models vs. the unspecialized parent model
    function relLikeli = relativeLikelihood(obj, gofFcn)
      %% Default arguments
      if ~exist('gofFcn', 'var') || isempty(gofFcn)
        gofFcn            = @CategorizedModel.poissonLogLikelihood;
      end
      
      %% Get the set of all unique parent models
      if any(~obj.hasProgenitor())
        error('CategorizedModel:relativeLikelihood', 'Only defined for model sets with parent models.');
      end
      parent              = unique([obj.progenitor]);
      
      %% Construct test statistic as the difference between model goodness of fit
      [setPredict, allTarget] = obj.jointPrediction();
      jointGOF            = gofFcn(allTarget, setPredict              , sum(obj.effectiveDegreesOfFreedom()));
      parentGOF           = gofFcn(allTarget, parent.jointPrediction(), sum(parent.effectiveDegreesOfFreedom()));
      
      %% Relative likelihood N.B. assuming that the goodness-of-fit is equivalent to the log-likelihood
      % A small value indicates that the specialized model is more probable than the parent model
      relLikeli           = exp(parentGOF - jointGOF);
    end
    
    %----- Sorts models in the given set, given a goodness-of-fit function
    function [obj, iOrder, gof] = sort(obj, fcn, direction, varargin)
      if nargin < 3
        direction               = '';
      elseif ~ischar(direction)
        varargin                = [{direction}, varargin];
        direction               = '';
      end
      
      if ischar(fcn)
        gof                     = obj.(fcn)(varargin{:});
        if ~isempty(direction)
        elseif numel(fcn) > 2 && strcmp(fcn(2:3), 'IC')
          direction             = 'ascend';
        else
          direction             = 'descend';
        end
      else
        gof                     = fcn(obj, varargin{:});
        if isempty(direction)
          direction             = 'ascend';
        end
      end
      
      [gof,iOrder]              = sort(gof, direction);
      obj                       = obj(iOrder);
    end
    
    
    %----- Regressor responses for each trial category
    function [response, time, numModels] = responseByCategory(obj, categoryLabels, categoryValues)
    % The Poisson GLM model predicts the neural responses in a given trial of category k to be
    % Poisson distributed with mean function
    %   \mu_k(t) = exp( m_k^{(1)}(t) + m_k^{(2)}(t) + ... + m_k^{(n)}(t) )
    % where the parenthesized indices i in m_k^{(i)} indicate temporal responses to the i-th type 
    % of behavioral event. If the specialization by trial category k is not of interest, i.e.
    % omitted from the categoryLabels input set, this function returns the average response
    %   \sum_k m_k^{(i)}(t)
    % for each behavioral event i, where the sum is over all N values k of the k-specialized models.
    % This is equivalent to taking the geometrical mean prediction \prod_k \mu_k(t) / N, as opposed
    % to the arithmetic mean.
      
      %% Check input format
      if ischar(categoryLabels)
        categoryLabels    = {categoryLabels};
      end
      if numel(categoryLabels) ~= size(categoryValues,2)
        error('CategorizedModel:responseByCategory', 'categoryValues must be an n-by-p matrix where n is the number of condition combinations to compute regressor responses for, and p is the number of categoryLabels.');
      end
      
      %% All models are required to be based on the same data
      experiment          = obj(1).design.dspec.expt;
      if any(~arrayfun(@(x) isequal(x.design.dspec.expt.id,experiment.id), obj(2:end)))
        error('CategorizedModel:responseByCategory', 'Inconsistent experimental data encountered for the given set of models.');
      end
      
      %% Preallocate output
      egRegressors        = obj(1).regressors;
      regType             = fieldnames(egRegressors);
      regType( cellfun(@(x) ~isstruct(egRegressors.(x)), regType) )  = [];
      numModels           = zeros(size(categoryValues,1), 1);
      response            = struct();
      time                = struct();
      for iReg = 1:numel(regType)
        response.(regType{iReg})  = zeros(numel(egRegressors.(regType{iReg}).response), size(categoryValues,1));
        time.(regType{iReg})      = egRegressors.(regType{iReg}).time;
      end
      
      %% Store category-specific regressors, taking the average across specializations for categories not in the desired categoryLabels
      specs               = obj.specializations();
      for iObj = 1:numel(obj)
        %% Identify categoryValues that this object contributes to
        isInCat           = true(size(categoryValues));
        for what = fieldnames(specs(iObj))'
          selCat          = strcmp(categoryLabels, what{:});
          isInCat(:,selCat) = ismember( categoryValues(:,selCat), specs(iObj).(what{:}) );
        end
        isInCat           = all(isInCat, 2);
        
        %% Add regressor contribution to all relevant categories
        for iReg = 1:numel(regType)
          for iCat = find(isInCat)'
            response.(regType{iReg})(:,iCat)      ...
                          = response.(regType{iReg})(:,iCat) + obj(iObj).regressors.(regType{iReg}).response;
            assert(isequal(time.(regType{iReg}), obj(iObj).regressors.(regType{iReg}).time));
          end
        end
        numModels         = numModels + isInCat;
      end
      
      %% Divide by number of contributing models to obtain average
      for iCat = 1:numel(response)
        if ~numModels(iCat)
          continue        % this can happen if the user provided categoryValues that are a superset of what the models were trained on
        end
        for iReg = 1:numel(regType)
          response(iCat).(regType{iReg})      = response(iCat).(regType{iReg}) / numModels(iCat);
        end
      end
    end
    
  end
  
end
