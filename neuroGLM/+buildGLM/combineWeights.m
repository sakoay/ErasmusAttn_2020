function [wout] = combineWeights(dm, w)
% Combine the weights per column in the design matrix per covariate
%
% Input
%   dm: design matrix structure
%   w: weight on the basis functions
%
% Output
%   wout.(label).basisW = weights for each basis function
%   wout.(label).response = combined weights
%   wout.(label).time = time axis

wout = struct();    %% SAK : prevent overwriting
dspec = dm.dspec;
binSize = dspec.expt.binSize;

if isfield(dm, 'biasCol') % undo z-score operation
    %% SAK : fixed biasCol addressing which must take into account the removed constCols
    biasCol = false(size(dm.constCols));
    biasCol(dm.biasCol) = true;
    biasCol(dm.constCols) = [];
    wout.bias = w(biasCol);
    w(biasCol) = [];
end

if isfield(dm, 'zscore') % undo z-score operation
    w = (w .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
end

if isfield(dm, 'constCols') % put back the constant columns
    w2 = zeros(dm.dspec.edim, 1);
    
    %% SAK : have to properly omit bias column
    constCols = dm.constCols;
    if isfield(dm, 'biasCol')
      constCols(dm.biasCol) = [];
    end
    
    w2(~constCols) = w; % first term is bias
    w = w2;
end

if numel(w) ~= dm.dspec.edim
    error('Expecting w to be %d dimension but it''s [%d]', ...
	dspec.edim, numel(w));
end

startIdx = [1 (cumsum([dspec.covar(:).edim]) + 1)];

for kCov = 1:numel(dspec.covar)
    covar = dspec.covar(kCov);
    basis = covar.basis;

    if isempty(basis)
      w_sub = w(startIdx(kCov) + (1:covar.edim) - 1);
      wout.(covar.label).basisW = w_sub;       % SAK
      wout.(covar.label).time = ((1:size(w_sub, 1))-1 + covar.offset) * binSize;
      wout.(covar.label).response = w_sub;
    	continue;
    end

    assert(isstruct(basis), 'Basis structure is not a structure?');

    sdim = covar.edim / basis.edim;
    wout.(covar.label).response = zeros(size(basis.B, 1), sdim);
    for sIdx = 1:sdim
      w_sub = w(startIdx(kCov) + (1:basis.edim)-1 + basis.edim * (sIdx - 1));
      w2_sub = sum(bsxfun(@times, basis.B, w_sub(:)'), 2);
      wout.(covar.label).basisW = w_sub;       % SAK
      wout.(covar.label).response(:, sIdx) = w2_sub;
    end
    wout.(covar.label).time = ...
	(basis.tr(:, 1) + covar.offset) * binSize * ones(1, sdim);
end
