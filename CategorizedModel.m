% Regression model for a subset of trials specified by some category selection criteria
classdef CategorizedModel < handle
  
  %_________________________________________________________________________________________________
  properties (SetAccess = protected)
    ID                = []            % unique identifier for each model, not changed upon copying
    
    progenitor        = []            % the parent model that includes data across all values of this model's category
    categoryLabel     = []            % the trial category label that is used to select data used in this model fit
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
    
    %----- Recursively find the maximum model ID value
    function id = maxID(obj)
      id      = max([obj.ID]);
      for iObj = 1:numel(obj)
        if obj(iObj).hasProgenitor()
          id  = max(id, obj(iObj).progenitor.maxID());
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
      if ~exist('prefix', 'var') || isempty(prefix)
        prefix        = ' ';
      end

      %% Print each model in this set
      for iObj = 1:numel(obj)
        fprintf('%s[%3d]  %*s%4d trials : ', prefix, obj(iObj).ID, 15-numel(prefix), '', obj(iObj).numberOfTrials());
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
          obj(iObj).progenitor.print([repmat(' ',size(prefix)), '\--']);
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
    
    %----- Get the set of all trial categories that this model is specialized for
    function labels = categories(obj, singletonOutput)
      labels            = cell(size(obj));
      for iObj = 1:numel(obj)
        if obj(iObj).hasProgenitor()
          labels{iObj}  = [obj(iObj).progenitor.categories(), {obj(iObj).categoryLabel}];
        elseif isempty(obj(iObj).categoryLabel)
          labels{iObj}  = {};
        else
          labels{iObj}  = {obj(iObj).categoryLabel};
        end
      end
      
      %% Special case for single objects
      if numel(obj) == 1 && (nargin < 2 || ~singletonOutput)
        labels          = labels{1};
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
        if ~isempty(obj(iObj).fitInfo)
          dof(iObj) = sum(obj(iObj).fitInfo.activeSet(:,obj(iObj).fitInfo.IndexMinDeviance), 1);
        end
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
    function obj = regress(obj, nCVFolds, nLambdas, lazy)
      %% Default arguments
      if ~exist('nLambdas', 'var') || isempty(nLambdas)
        nLambdas                      = 10;
      end
      if ~exist('lazy', 'var') || isempty(lazy)
        lazy                          = false;
      end

      %% GLM configuration
      fitOpts                         = {'poisson', 'Link', 'log', 'MaxIter', 1e5};
      
      %% Fit model for each configuration of model specifications and data to explain
      for iObj = 1:numel(obj)
        if lazy && ~isempty(obj(iObj).prediction)
          continue;
        end
        
        %% Define cross-validation selection and shuffled data
        trialLength                   = SplitVec(obj(iObj).design.trialX, 'equal', 'length');
        [~,~,obj(iObj).cvTrainSel]    = permutationExperiments(trialLength, 0, 0, nCVFolds, 2);
        xvalPartitions                = cvpartition(numel(trialLength), 'KFold', nCVFolds);
                                
        designX                       = full(obj(iObj).design.X);
        targetY                       = obj(iObj).target;
%         fitglm(designX, targetY, 'linear', 'Distribution', 'poisson', 'Link', 'log', 'Intercept', true );
        
%         fitInfo                   = linearSupportVectorMR( obj(iObj).target, obj(iObj).design.X, obj(iObj).cvTrainSel, obj(iObj).shuffleExp, bootstrapExp, [], [], [], true );

        %% Use cross-validation to determine regularization strength
        xvalPartitions.Impl           = CustomCVPartition(obj(iObj).cvTrainSel(:,:,1));
        [modelW,obj(iObj).fitInfo]    = lassoglm(designX, targetY, fitOpts{:}, 'NumLambda', nLambdas, 'CV', xvalPartitions);

        lambdaRatio                   = obj(iObj).fitInfo.LambdaMinDeviance / obj(iObj).fitInfo.Lambda(end);
%         obj(iObj).fitInfo.LambdaL1   = obj(iObj).fitInfo.Lambda .* obj(iObj).fitInfo.Alpha;
%         obj(iObj).fitInfo.LambdaL2   = obj(iObj).fitInfo.Lambda .* (1 - obj(iObj).fitInfo.Alpha)/2;
        obj(iObj).fitInfo.activeSet   = modelW ~= 0;
        
%         lassoPlot(cvW,cvFitInfo,'plottype','CV');

        %% Use an independent set of cross-validation folds to evaluate goodness of fit
        if lambdaRatio < 1
          cvPrediction                = nan(size(obj(iObj).target));
          for iCV = 1:nCVFolds
            %% Train model on a subset of data
            trainSel                  = obj(iObj).cvTrainSel(:,iCV,2);
            [cvW,cvInfo]              = lassoglm(designX, targetY, fitOpts{:}, 'NumLambda', 2, 'LambdaRatio', lambdaRatio);

            %% Evaluate model on left-out test data
            testSel                   = ~trainSel;
            cvPrediction(testSel)     = glmval([cvInfo.Intercept(1); cvW(:,1)], designX(testSel,:), 'log');
          end
        else
          cvPrediction                = zeros(size(obj(iObj).target));
        end
        
        %% Store fitted parameters
        obj(iObj).prediction          = cvPrediction;
        obj(iObj).regressors          = buildGLM.combineWeights(obj(iObj).design, modelW(:,obj(iObj).fitInfo.IndexMinDeviance));
        obj(iObj).regressors.bias     = obj(iObj).fitInfo.Intercept(obj(iObj).fitInfo.IndexMinDeviance);
%         obj(iObj).shufflePvalue   = mean( deviance <= deviance(1) );    % N.B. the unshuffled experiment is included as a pseudo-count to provide a conservative nonzero estimate
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
          gof(:,iObj)     = gofFcn(obj(iObj).target, obj(iObj).prediction, obj(iObj).fitInfo.activeSet);
        else
          gof(iObj)       = gofFcn(obj(iObj).target, obj(iObj).prediction, obj(iObj).fitInfo.activeSet);
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
    
    %----- Assign rankings (1 is best) to all models in the given set, given a goodness-of-fit function 
    function [obj, gof, iOrder, ranking] = assignRanks(obj, fcn, direction, varargin)
      if nargin < 3
        direction               = '';
      elseif ~isempty(direction) && ~ischar(direction)
        varargin                = [{direction}, varargin];
        direction               = '';
      end
      [obj, iOrder, gof]        = obj.sort(fcn, direction, varargin{:});
      ranking                   = num2cell(1:numel(obj));
      [obj.ranking]             = ranking{:};
    end
  end
  
end
