function [L,corrval,overlap] = vloss(XTest, R, M, V, delta,reg,hypers)
% normal loss function based on condition-specific variances
%
% XTest  testing dataset
% R      correlation matrix of training set
% Rp     correlation matrix of testing set based on training variances
% V      nBins*nConds*1*nCells - variances per bin per cell
% N      nBins*nConds*1*nCells - numbers of trials for each element of V in XTest



p = size(XTest,4);
assert(size(V,4)==p)
N = sum(~isnan(XTest),3);
N = N/sum(N(:));
assert(all(size(N)==size(V)))

% subtract means
XTest = bsxfun(@minus, XTest, M);

% shrink the variance estimates toward the mean variance across all conditions
if delta
    V0 = reshape(nanmean(reshape(V,[],p),1),[1 1 1 p]);
    V = bsxfun(@plus,(1-delta)*V,delta*V0);
end

% normalize by regularized variance
Z = bsxfun(@rdivide, XTest, sqrt(V));

% covariance of z-score, i.e. the correlation estimate
Rp = cove.cov(reshape(Z,[],p));

% normal loss
L = (trace(Rp/R)+cove.logDet(R))/p + sum(log(V(:)).*N(:));

corrval = [];
overlap = [];

if nargin>5
    % correlation between test and train correlation matrices
    corrval = corr(squareform(nodiag(Rp))',squareform(nodiag(R))');
    
    overlap = edgeoverlap(R,Rp,5);
    
%     sigma = diag(sqrt(mean(reshape(V,[],p))));  % average variances
%     C = (sigma*Rp*sigma);
%     C = cove.estimate_simple(C,reg, hypers); % get model-based covariance
%     Rp_reg = sigma\C/sigma;
        
end

if isnan(L)
    L = inf;
end