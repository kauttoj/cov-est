function [DATA,corrmat,covmat] = make_random_network(nCells,nTrials,Bins,nConditions)
%MAKE_RANDOM_NETWORK Summary of this function goes here
%   Detailed explanation goes here

%%
if nargin==0
    
    nCells=50;
    nTrials=1000;
    Bins=5;
    nConditions = 1;
    
end

nCommon = 3;
covmat = 0;

for k=1:nConditions
    
    A = 2*rand(nCommon,nCells).^4;
    common_signals = 5*randn(nCommon,Bins)';
    
    SIGMA = triu(randn(nCells,nCells),1);
    SIGMA(abs(SIGMA)<2.5)=0;
    SIGMA = SIGMA + SIGMA'  + 1.5*rand()*diag(rand(1,nCells));
    SIGMA = SIGMA*SIGMA';
        
    signals = nan(Bins,nCells,nTrials);
    for i=1:nTrials
        signals(:,:,i) = common_signals*A + mvnrnd(zeros(1,nCells),SIGMA,Bins);
    end
    
    DATA(:,:,:,k) = signals;    
    
    covmat = covmat + SIGMA/nConditions;
    
end

%X - data: nBins * nConds * nTrials * nCells

corrmat = corrcov(SIGMA);

DATA = permute(DATA,[1,4,3,2]);

end

