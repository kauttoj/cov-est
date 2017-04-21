% simulation testbed

clearvars
close all;
clc;

rng(666);
% nCells = 200;
% T = 200; % trials/repeats
% Bins = 2;%
% 
% %%
% 
% t = linspace(0,50*pi,T*Bins);
% COMMON1 = sin(0.7*t+0.5).*(cos(t).^2) + 0.05*randn(1,T*Bins);
% COMMON2 = sin(0.3*sqrt(t)+0.5).*(cos(1.1*t).^2) + 0.2*rand(1,T*Bins);
% 
% A = 0.2+rand(2,nCells);
% A(1,:)=A(1,:)/2;
% 
% SIGMA = triu(randn(nCells,nCells),1);
% SIGMA(SIGMA<2.5)=0;
% SIGMA = 3*SIGMA;
% SIGMA = SIGMA + SIGMA'  + 3*diag(rand(1,nCells));
% SIGMA = SIGMA*SIGMA';
% 
% covmat_true = SIGMA;
% cormat_true = corrcov(SIGMA);
% 
% signals = mvnrnd(zeros(1,nCells),SIGMA,T*Bins);
% 
% for i=1:nCells
%     signals(:,i) = signals(:,i) + 10*A(1,i)*COMMON1'+ 5*A(2,i)*COMMON2';
% end
% 
% DATA = signals;
% DATA = reshape(DATA,[Bins,T,nCells]);
% DATA = permute(DATA,[1,4,2,3]);
% 
% DATA_small = DATA(:,:,round(size(DATA,3)/2),:);
% 
% IND = find(triu(ones(nCells,nCells),1));


%X - data: nBins * nConds * nTrials * nCells

SAMPLES = 100;
CELLS = 200;
BINS = 5;
CONDITIONS = 2;
GRID = {logspace(-4,-0.2,12),logspace(-4,-0.2,12)};

[DATA,real_corrmat,real_covmat] = make_random_network(CELLS,SAMPLES,BINS,CONDITIONS);

%% REAL DATA
% A = load('traces-pre5-0001.mat');
% DATA = A.tuples(23).trace_segments;
% DATA = DATA(1:BINS,[1,4,5],:,randsample(size(DATA,4),CELLS));
% real_corrmat = zeros(CELLS,CELLS);
% real_corrmat(1,1)=1;
%%%

IND = find(triu(ones(CELLS,CELLS),1));

figure('position',[100         892        2102         487]);

subplot(1,5,1);imagesc(nodiag(real_corrmat));axis image;colorbar;
title('real')

[hypers, bestDelta, visited, losses, quality1] = cove.crossEstimateHyper(DATA,size(DATA,1),'sample');%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat1,M1, V1, extras1,covmat1] = cove.estimate(DATA,size(DATA,1),'sample',[]); % alpha = 0.05; beta = 0.25
subplot(1,5,2);imagesc(nodiag(cormat1));axis image;colorbar;
title(sprintf('sample (%f)',corr(real_corrmat(IND),cormat1(IND))));

[hypers, bestDelta, visited, losses, quality4] = cove.crossEstimateHyper(DATA,size(DATA,1),'diag',GRID);%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat4,M4, V4, extras4,covmat4] = cove.estimate(DATA,size(DATA,1),'diag',hypers); % alpha = 0.05; beta = 0.25
subplot(1,5,5);imagesc(nodiag(cormat4));axis image;colorbar;
title(sprintf('yatsenko diag (%f)',corr(real_corrmat(IND),cormat4(IND))));

[hypers, bestDelta, visited, losses, quality3] = cove.crossEstimateHyper(DATA,size(DATA,1),'lv-glasso',GRID);%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat3,M3, V3, extras3,covmat3] = cove.estimate(DATA,size(DATA,1),'lv-glasso',hypers); % alpha = 0.05; beta = 0.25
subplot(1,5,4);imagesc(nodiag(cormat3));axis image;colorbar;
title(sprintf('yatsenko model (%f)',corr(real_corrmat(IND),cormat3(IND))));

[hypers, bestDelta, visited, losses, quality2] = cove.crossEstimateHyper(DATA,size(DATA,1),'glasso',GRID(1));%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat2,M2, V2, extras2,covmat2] = cove.estimate(DATA,size(DATA,1),'glasso',hypers); % alpha = 0.05; beta = 0.25
subplot(1,5,3);imagesc(nodiag(cormat2));axis image;colorbar;
title(sprintf('yatsenko glasso (%f)',corr(real_corrmat(IND),cormat2(IND))));


figure('position',[100        100        2102         487]);
subplot(1,4,1);imagesc(nodiag(real_corrmat));axis image;colorbar;
title('real')

a = -corrcov(inv(covmat1));
subplot(1,4,2);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('sample sparse (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));

a = -corrcov(extras2.S);
subplot(1,4,3);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('yatsenko model sparse (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));

a = -corrcov(inv(covmat3));
subplot(1,4,4);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('yatsenko diag sparse (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));



