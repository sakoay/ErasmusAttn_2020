function [dm, trialBins, eventEpochs] = compileSparseDesignMatrix(dspec, trialIndices)
% Compile information from experiment according to given DesignSpec

%% Build design matrix
expt = dspec.expt;
subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);

totalT = sum(ceil([expt.trial(trialIndices).duration]/expt.binSize));

growingX = sparse([], [], [], 0, dspec.edim, round(totalT * dspec.edim * 0.001)); % preallocate
trialX = nan(totalT,1);                 % SAK : store index of the trial corresponding to a given row of the design matrix
trialBins = cell(1,max(trialIndices));  % SAK : store rows of the design matrix corresponding to the input set of trials
eventEpochs = nan(totalT,1);            % SAK : store linearized time (epochs) between events in each trial
% lastEpoch = linspace(numel(dspec.covar), numel(dspec.covar)+1, size(dspec.covar(end).basis.B,1) + 1);

trialIndices = trialIndices(:)';

for kTrial = trialIndices
    nT = ceil(expt.trial(kTrial).duration / expt.binSize); % TODO move?
    
    miniX = zeros(nT, dspec.edim);  % pre-allocate a dense matrix for each trial
    miniEpoch = nan(nT, 1);         % SAK
    prevCov = [];                   % SAK
    prevIdx = [];                   % SAK
    
    for kCov = 1:numel(dspec.covar) % for each covariate
        covar = dspec.covar(kCov);
        sidx = subIdxs{kCov};
        
        if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(kTrial))
            continue;
        end
        
        stim = covar.stim(expt.trial(kTrial), nT); % either dense or sparse
        
        if isfield(covar, 'basis') && ~isempty(covar.basis)
            miniX(:, sidx) = basisFactory.convBasis(stim, covar.basis, covar.offset);
            
            %% SAK
            startIdx = find(stim, 1, 'first');
            if ~isempty(startIdx)
              if ~isempty(prevCov)
                miniEpoch(prevIdx:startIdx) = linspace(prevCov, kCov, startIdx - prevIdx + 1);
              end
              prevCov = kCov;
              prevIdx = startIdx;
            end
        else
            miniX(:, sidx) = stim;
        end
    end
    
    trialBins{kTrial} = size(growingX,1) + (1:nT);  % SAK
    trialX(trialBins{kTrial}) = kTrial;             % SAK
    growingX = [growingX; sparse(miniX)]; %#ok<AGROW>
    
    %% SAK -- UGLY : keep linear time past the last event
    if ~isempty(prevIdx)
      lastEpoch = linspace(prevCov, numel(dspec.covar)+1, (numel(dspec.covar)+1 - prevCov)*size(dspec.covar(end).basis.B,1) + 1);
      epochIdx = 1:min(size(miniEpoch,1) - prevIdx + 1, numel(lastEpoch));
      miniEpoch(prevIdx-1 + epochIdx) = lastEpoch(epochIdx);
    end
    eventEpochs(trialBins{kTrial}) = miniEpoch;
end

dm.X = growingX;
dm.trialX = trialX;
dm.trialIndices = trialIndices;
dm.dspec = dspec;

%% Check sanity of the design
if any(~isfinite(dm.X(:)))
    warning('Design matrix contains NaN or Inf...this is not good!');
end
