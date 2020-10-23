function dm = compileSparseDesignMatrix(dspec, trialIndices, trialCategory)
% Compile information from experiment according to given DesignSpec

%% SAK : support duplication of covariates so that there is a different set associated with each indexed trial category
if nargin < 3 || isempty(trialCategory)
  categIndex = [];
  nCategories = 1;
else
  if size(trialCategory,1) < max(trialIndices)
    error('compileSparseDesignMatrix:trialCategory', 'trialCategory must have rows that can be indexed by trialIndices.');
  end
  
  %% Identify the unique category to which the specified trials belong
  [category, ~, categIndex] = unique(trialCategory, 'rows');
  nCategories = size(category,1);
end

%% Build design matrix
expt = dspec.expt;
subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);

totalT = sum(ceil([expt.trial(trialIndices).duration]/expt.binSize));

growingX = sparse([], [], [], 0, nCategories*dspec.edim, round(totalT * dspec.edim * 0.001)); % preallocate
trialX = nan(totalT,1);   % SAK : store index of the trial corresponding to a given row of the design matrix

trialIndices = trialIndices(:)';

for kTrial = trialIndices
    nT = ceil(expt.trial(kTrial).duration / expt.binSize); % TODO move?
    
    miniX = zeros(nT, nCategories*dspec.edim); % pre-allocate a dense matrix for each trial
    
    for kCov = 1:numel(dspec.covar) % for each covariate
        covar = dspec.covar(kCov);
        sidx = subIdxs{kCov};
        
        if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(kTrial))
            continue;
        end
        
        if ~isempty(categIndex)
          sidx = sidx + dspec.edim * (categIndex(kTrial) - 1);    % SAK covariates by trial category
        end
        stim = covar.stim(expt.trial(kTrial), nT); % either dense or sparse
        
        if isfield(covar, 'basis') && ~isempty(covar.basis)
            miniX(:, sidx) = basisFactory.convBasis(stim, covar.basis, covar.offset);
        else
            miniX(:, sidx) = stim;
        end
    end
    trialX(size(growingX,1) + (1:nT)) = kTrial;   % SAK
    growingX = [growingX; sparse(miniX)]; %#ok<AGROW>
end

dm.X = growingX;
dm.trialX = trialX;
dm.trialIndices = trialIndices;
if ~isempty(categIndex)
  dm.categories = category;
  dm.trialCategory = categIndex(trialIndices);
  dm.trialCategory = dm.trialCategory(:)';
end
dm.dspec = dspec;

%% Check sanity of the design
if any(~isfinite(dm.X(:)))
    warning('Design matrix contains NaN or Inf...this is not good!');
end
