function bases = makeSmoothTemporalBasis(shape, duration, nBases, binfun)
%
% Input
%   shape: 'raised cosine' or 'boxcar'
%   duration: the time that needs to be covered
%   nBases: number of basis vectors to fill the duration
%   binfun: 
%
% Output
%   BBstm: basis vectors

nkbins = binfun(duration); % number of bins for the basis functions

ttb = repmat((1:nkbins)', 1, nBases); % time indices for basis

if strcmpi(shape, 'raised cosine')
    %   ^
    %  / \
    % /   \______
    %      ^
    %     / \
    % ___/   \___
    %         ^
    %        / \
    % ______/   \
    % For raised cosine, the spacing between the centers must be 1/4 of the
    % width of the cosine
    dbcenter = nkbins / (3 + nBases); % spacing between bumps
    width = 4 * dbcenter; % width of each bump
    bcenters = 2 * dbcenter + dbcenter*(0:nBases-1);
    % location of each bump centers
    % bcenters = round(bcenters); % round to closest ms <-- BAD!!! NEVER DO THIS
    bfun = @(x,period)((abs(x/period)<0.5).*(cos(x*2*pi/period)*.5+.5));
    BBstm = bfun(ttb-repmat(bcenters,nkbins,1), width);
elseif strcmpi(shape, 'progressive cosine x2')
    %% SAK : cosine bumps with progressively x2 increases in width and spacing, with the first component at zero
    sigma = nkbins * 2 / (2^nBases - 1); % width of first bump
    width = 2.^(0:nBases-1) * sigma; % width of each bump
    bcenters = 1 + [0, (2.^(2:nBases) - 2)*sigma/4];    % location of each bump center
    bfun = @(x,period)((abs(x./period)<0.5).*(cos(x*2*pi./period)*.5+.5));
    BBstm = bfun(ttb-bcenters, width);
%     figure; plot(BBstm)
elseif strcmpi(shape, 'progressive cosine x1.5')
    %% SAK : cosine bumps with progressively x1.5 increases in width and spacing, with the first component at zero
    sigma = nkbins / ((3/2)^nBases - 1); % width of first bump
    width = (3/2).^(0:nBases-1) * sigma; % width of each bump
    bcenters = 1 + [0, ((3/2).^(1:nBases-1) - 1)*sigma];    % location of each bump center
    bfun = @(x,period)((abs(x./period)<0.5).*(cos(x*2*pi./period)*.5+.5));
    BBstm = bfun(ttb-bcenters, width);
    
%     longfigure(2,0.5); plot((0:nkbins-1) * duration / nkbins, BBstm,'linewidth',1.5); xlabel('Time from event (ms)'); ylabel('Basis function')

elseif strcmpi(shape, 'progressive cosine x2 off0')
    %% SAK : cosine bumps with progressively x2 increases in width and spacing, with the first component offset from zero
    sigma = nkbins * 4 / (3*2^nBases - 2); % width of first bump
    width = 2.^(0:nBases-1) * sigma; % width of each bump
    bcenters = 1 + (2.^(1:nBases) - 1)*sigma/2;    % location of each bump center
    bfun = @(x,period)((abs(x./period)<0.5).*(cos(x*2*pi./period)*.5+.5));
    BBstm = bfun(ttb-bcenters, width);
    
elseif strcmpi(shape, 'progressive cosine x1.5 off0')
    %% SAK : cosine bumps with progressively x1.5 increases in width and spacing, with the first component offset from zero
    sigma = nkbins * 3 * 2^nBases / (4*3^nBases - 3*2^nBases); % width of first bump
    width = (3/2).^(0:nBases-1) * sigma; % width of each bump
    bcenters = 1 + ((3/2).^(1:nBases) - 1)*sigma;    % location of each bump center
    bfun = @(x,period)((abs(x./period)<0.5).*(cos(x*2*pi./period)*.5+.5));
    BBstm = bfun(ttb-bcenters, width);
    
elseif strcmpi(shape, 'boxcar')
    width = nkbins / nBases;
    BBstm = zeros(size(ttb));
    bcenters = width * (1:nBases) - width/2;
    for k = 1:nBases
        idx = ttb(:, k) > ceil(width * (k-1)) & ttb(:, k) <= ceil(width * k);
        BBstm(idx, k) = 1 / sum(idx);
    end
else
    error('Unknown basis shape');
end

bases.type = [shape '@' mfilename];
bases.param.shape = shape;
bases.param.duration = duration;
bases.param.nBases = nBases;
bases.param.binfun = binfun;
bases.B = BBstm;
bases.edim = size(bases.B, 2);
bases.tr = ttb - 1;
bases.centers = bcenters;
