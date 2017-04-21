function [R, M, V, extras,C] = estimate(X, evokedBins, reg, hypers, PRINT_STATS)
% estimate covariance matrix C
%
% Input:
%    X - data: nBins * nConds * nTrials * nCells
%    evokedBins - the number of bins considered in the evoked response,
%    whereas the remaning bins are considered to be spontaneous
%    reg - covariance regularization method
%
% Output:
%    R  - correlation matrix
%    M  - mean responses:  nBins * nConds * nTrials * nCells
%    V  - variances, dimensions:  nBins * nConds * 1 * nCells
%    extras - structure with additional information about the estimate,
%    i.e. sparse component, low-rank components, sparsity, etc.

warning('off','MATLAB:nearlySingularMatrix');

if nargin==0
   X = randn(7, 3, 300,50);
end

if nargin<5
   PRINT_STATS = 1; 
end

X = double(X);

extras = struct;
[nBins, nConds, nTrials, nCells] = size(X);
evokedBins = nBins;

% estimate mean response
% bin-specific means
M = nanmean(X(1:evokedBins,:,:,:),3);   % binwise mean of evoked response

% subtract mean
X = bsxfun(@minus, X, M);

% compute variance over bins instead of separately for bins
XX = permute(X,[1,3,2,4]);
XX = reshape(XX,[nBins*nTrials,nConds,nCells]);
V = nanvar(XX,1,1);
V = permute(V,[1,2,4,3]);
V = repmat(V,[nBins,1,1,1]);

% estimate binwise variances
%V = nanvar(X(1:evokedBins,:,:,:),1,3);   % binwise mean of evoked response

% sample correlation matrix based on bin-wise variances
Z = bsxfun(@rdivide, X, sqrt(V));
R = cove.cov(reshape(Z,[],nCells));
assert(all(abs(diag(R)-1)<1e-3))
R = corrcov(R,true);  % just in case

% average variances across all bins and produce the average
% sample covariance matrix, C
sigma = diag(sqrt(mean(reshape(V,[],nCells))));  % average variances
C = (sigma*R*sigma);

R_sample = R;

switch reg
    case 'sample'
        assert(isempty(hypers),'invalid hyperparameters')
        % do nothing
        
    case 'diag'
        % 2 hypers:  variance shrinkage,  correlation shrinkage
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        r = (1-hypers(2))*corrcov(C) + hypers(2)*eye(size(C));
        C = v*r*v;
        
    case 'factor'
        % 2 hypers: variance shrink toward median, nlatent
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        C = v*corrcov(C)*v;
        [L,psi] = cove.factor(C,hypers(2));  % factor analysis
        extras.loading_matrix = L;
        extras.indep_vars = psi(:);
        C = L*L' + diag(psi);
        
    case 'glasso'
        assert(length(hypers)==1)
        cove.set('max_latent',0)   % prevent latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),10,cove.set);
        extras.S = extras.S/scale;
        C = inv(extras.S);
        
    case 'L1'        
        assert(length(hypers)==1)
        C = L1precisionBCD(C/mean(diag(C)),hypers(1));
        
    case 'lv-glasso'
        assert(length(hypers)==2)
        cove.set('max_latent',inf)   % allow any number of latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),hypers(2),cove.set);
        extras.S = extras.S/scale;  % scale back
        extras.L = extras.L/scale;
        [H,D] = svds(double(extras.L),sum(~~extras.eigL));
        extras.H = H*sqrt(D);
        C = inv(extras.S - extras.H*extras.H');        

    otherwise
        error 'unknown covariance estimator'
end

% convert back to correlations
R = sigma\C/sigma;

C = (C+C')/2;

extras.CondEst = condest(C);

if ~strcmp(reg,'sample')
    % transfer the change in variance from R to V
    V = bsxfun(@times, V, reshape(diag(R), [1 1 1 nCells]));
    R = corrcov(R);    
end

if PRINT_STATS
    if ~strcmp(reg,'sample')
        fprintf('corr2(R_sample,R_reg) = %f, condest(C_reg) = %f, type ''%s''\n\n',corr2(squareform(nodiag(R_sample))',squareform(nodiag(R))'),extras.CondEst,reg);
    else
        fprintf('condest(C_sample) = %f, type ''%s''\n\n',extras.CondEst,reg);
    end
end