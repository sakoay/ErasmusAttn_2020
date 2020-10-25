% Regression model for a subset of trials specified by some category selection criteria
classdef CategorizedModel
  
  %_________________________________________________________________________________________________
  properties (SetAccess = protected)
    ID                = []            % unique identifier for each model, not changed upon copying
    
    progenitor        = []            % the parent model that includes data across all values of this model's category
    categoryLabel     = []            % the trial category label that is used to select data used in this model fit
    categoryValues    = []            % values of categoryLabel that are included in the selected data in the fit for this model

    target            = []            % the time series data predicted by this model
    design            = []            % design matrix specifying covariates used to predict target 
    prediction        = []            % cross-validated prediction for target 
    model             = []            % fitted model
    regressors        = struct()      % fitted weights, time axis and waveforms for regressors

    shuffleExp        = []            % shuffled time indices used to generate null hypotheses where there is no relationship between target and regressors
    cvTrainSel        = []            % selection vector for subsets of target used to generate training sets for cross-validation
    altModels         = {}            % alternative models that were considered
    extra             = struct()      % optional user specified information
  end
  
  %_________________________________________________________________________________________________
  methods (Static)
    
    %----- Apply additional trial selections to design matrix
    function [design, selData] = selectTrials(design, selTrials)
      design.trialIndices(~selTrials)     = [];

      % Select timepoints corresponding to the selected trials
      selData                             = ismember(design.trialX, design.trialIndices);
      design.X(~selData,:)                = [];
      design.trialX(~selData,:)           = [];

      % Remove constant columns of the new design matrix
      isConst                             = full(all(design.X(1,:) == design.X(2:end,:), 1));
      design.constCols(~design.constCols) = isConst;
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
      if nargin > 0 && reset
        tally   = 0;
      elseif isempty(tally)
        tally   = 1;
      else
        tally   = tally + 1;
      end
      id        = tally;
    end

    %----- Variance explained: R^2 = 1 - SS_res / SS_tot where SS_x is the sum-squared values
    function R2 = computeVarExplained(target, prediction)
      sel     = isfinite(target);
      SSres   = sum( (target(sel) - prediction(sel)  ).^2, 1 );
      SStot   = sum( (target(sel) - mean(target(sel))).^2, 1 );
      R2      = 1 - SSres / SStot;
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

    
    %----- Create more specialized models from all possible partitions of the data by values of a given category
    function subModels = specialize(obj, categoryLabel)

      if numel(obj) > 1
        error('CategorizedModel:specialize', 'This function should be called on individual CategorizedModel objects.');
      end
      
      %% Identify all partitions of values of dataCategory in the dataset for this model
      dataCategory                = cat(1, obj.design.dspec.expt.trial.(categoryLabel));
      category                    = dataCategory(obj.design.trialIndices);
      categValues                 = unique(category);
      categorySets                = partitions(categValues);
      categorySets( cellfun(@numel, categorySets) < 2 ) = [];     % disallow trivial partitions
        
      %% Create models per partition
      subModels                   = repmat({repmat(CategorizedModel,0)}, size(categorySets));
      for iSet = 1:numel(categorySets)
        for iVal = 1:numel(categorySets{iSet})
          %% Apply further trial selection to design matrix
          [design, selData]       = CategorizedModel.selectTrials(obj.design, ismember(category, categorySets{iSet}{iVal}));
          
          %% Construct a more specialized model
          model                   = CategorizedModel(obj.target(selData), design);
          model.progenitor        = obj;
          model.categoryLabel     = categoryLabel;
          model.categoryValues    = categorySets{iSet}{iVal};
          subModels{iSet}(end+1)  = model;
        end
      end
      
    end
    
    %----- Fit regression model
    function obj = regress(obj, nShuffles, nCVFolds)
      %% Fit model for each configuration of model specifications and data to explain
      for iObj = 1:numel(obj)
        %% Define cross-validation selection and shuffled data
        trialLength               = SplitVec(obj(iObj).design.trialX, 'equal', 'length');
        [obj(iObj).shuffleExp, bootstrapExp, obj(iObj).cvTrainSel]       ...
                                  = permutationExperiments(trialLength, nShuffles, 0, nCVFolds);
        data                      = cell2table( num2cell([obj(iObj).design.X, obj(iObj).target])      ...
                                              , 'VariableNames', [accumfun(2, @(x) arrayfun(@(y) sprintf('%s_b%d',x.label,y), 1:x.edim, 'uniformoutput', false), obj(iObj).design.dspec.covar), {'firing_rate'}] );
                                            
%         fitInfo                   = linearSupportVectorMR( obj(iObj).target, obj(iObj).design.X, obj(iObj).cvTrainSel, obj(iObj).shuffleExp, bootstrapExp, [], [], [], true );
        
        %% Stepwise regression model to select a parsimonious subset of covariates to explain data
        obj(iObj).model           = stepwiseglm( data, 'constant', 'Criterion', 'Deviance', 'Distribution', 'poisson', 'Link', 'log', 'Intercept', true, 'Upper', 'linear', 'Verbose', 0 );
        
        %% Retrieve coefficient matrix with zeros for non-selected regressors
        [~,varIndex]              = ismember(obj(iObj).model.CoefficientNames(2:end), data.Properties.VariableNames);
        modelW                    = zeros(size(data,2), 1);
        modelW(1)                 = obj(iObj).model.Coefficients{'(Intercept)', 'Estimate'};
        modelW(varIndex)          = obj(iObj).model.Coefficients{2:end, 'Estimate'};
        
        %% Cross-validated predictions
%         cvW                       = nan(size(obj(iObj).design.X,2) + 1, size(obj(iObj).shuffleExp,2));
        cvPrediction              = nan(size(obj(iObj).target));
        for iCV = 1:size(obj(iObj).cvTrainSel,2)
          %% Train model on a subset of data
          trainSel                = obj(iObj).cvTrainSel(:,iCV);
%           [cvW(:,iCV),~,stats]    = glmfit(obj(iObj).design.X(trainSel,:), obj(iObj).target(trainSel), 'poisson', 'link', 'log');
          trainModel              = fitglm( data(trainSel,[varIndex, end]), 'linear', 'Distribution', 'poisson', 'Link', 'log', 'Intercept', true );
          
          %% Evaluate model on left-out test data
          testSel                 = ~trainSel;
%           cvPrediction(testSel)   = glmval(cvW(:,iCV), obj(iObj).design.X(testSel,:), 'log');
          cvPrediction(testSel)   = trainModel.predict(data(testSel,varIndex));
        end
        
        %% Store fitted parameters
        obj(iObj).regressors      = buildGLM.combineWeights(obj(iObj).design, modelW(2:end,1));
        obj(iObj).regressors.offset = modelW(1,1);
        obj(iObj).prediction      = cvPrediction;
%         obj(iObj).shufflePvalue   = mean( deviance <= deviance(1) );    % N.B. the unshuffled experiment is included as a pseudo-count to provide a conservative nonzero estimate
      end
    end
    
    
    %----- Variance explained: R^2 = 1 - SS_res / SS_tot where SS_x is the sum-squared values
    function gof = varExplained(obj)
      %% Preallocate output 
      nExperiments      = max(arrayfun(@(x) size(x.prediction,2), obj));
      if nExperiments > 1
        gof             = nan(nExperiments, numel(obj));
      else
        gof             = nan(size(obj));
      end
      
      %% Loop over models
      for iObj = 1:numel(obj)
        if isempty(obj(iObj).prediction)
          continue;
        end
        
        %% Coefficient of variation
        R2              = CategorizedModel.computeVarExplained(obj(iObj).target, obj(iObj).prediction);
        if nExperiments > 1
          gof(:,iObj)   = R2;
        else
          gof(iObj)     = R2;
        end
      end
    end
    
    %----- Cohen's f^2 for effect size of model A vs. (reference) B: (R^2_A - R^2_B) / (1 - R^2_A)
    function effect = relativeEffectSize(obj, reference)
      if isnumeric(reference)
        refR2             = reference;
      else
        refR2             = reference.varExplained();
      end
      effect              = obj.varExplained();
      effect              = (effect - refR2) ./ (1 - effect);
    end
    
    %----- Cohen's f^2 for local effect size vs. progenitor: (R^2_model - R^2_progen) / (1 - R^2_model)
    % Selya, A. S., Rose, J. S., Dierker, L. C., Hedeker, D. & Mermelstein, R. J. 
    % A Practical Guide to Calculating Cohen’s f2, a Measure of Local Effect Size, from PROC MIXED.
    % Front. Psychol. 3, (2012).
    function effect = localEffectSize(obj)
      %% All models are required to have the same parent model
      parent              = [obj.progenitor];
      if any(~obj.hasProgenitor())
        error('CategorizedModel:localEffectSize', 'Only defined for models with a parent model.');
      elseif any(parent(1).ID ~= [parent(2:end).ID])
        error('CategorizedModel:localEffectSize', 'Inconsistent parent models encountered.');
      end
      parent              = parent(1);
      
      %% The union of all models in this set are required to explain the same data as their parent model
      trials              = accumfun(1, @(x) x.trialIndices(:), [obj.design]);
      if numel(trials) ~= numel(parent.design.trialIndices) || ~isempty(setdiff(parent.design.trialIndices, trials))
        error('CategorizedModel:localEffectSize', 'Union of models in the given set must explain the same data as the parent model.');
      end
      
      %% Collect predictions from all models
      jointPredict        = nan(size(parent.target));
      for iObj = 1:numel(obj)
        subData           = ismember(parent.design.trialX, obj(iObj).design.trialIndices);
        jointPredict(subData) = obj(iObj).prediction;
      end
      R2                  = CategorizedModel.computeVarExplained(parent.target, jointPredict);
      
      %% Cohen's f^2
      effect              = (R2 - parent.varExplained()) / (1 - R2);
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
