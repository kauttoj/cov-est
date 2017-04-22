% simulation testbed

clearvars
close all;
clc;

rng(666);

%X - data: nBins * nConds * nTrials * nCells

% SAMPLES = 100;
CELLS = 230;
BINS = 5;
CONDITIONS = 2;
% [DATA,real_corrmat,real_covmat] = make_random_network(CELLS,SAMPLES,BINS,CONDITIONS);

%% REAL DATA
% A = load('data\traces-pre5-0001.mat');
% DATA = A.tuples(23).trace_segments;
% DATA = DATA(1:BINS,[1,4,5],:,randsample(size(DATA,4),CELLS));

%A=load('data\2093_NC_170408_data.mat');
A=load('data\2093_NC_170414_data.mat');
DATA = permute(A.DATA,[3,2,4,1]);
DATA = DATA(:,1:CONDITIONS,:,randsample(size(DATA,4),min(size(DATA,4),CELLS)));
%DATA = DATA + (1e-6)*rand(size(DATA)); % just in case there are null responses in one fold

real_corrmat = zeros(CELLS,CELLS);
real_corrmat(1,1)=1;

%%%

%alpha = 0.002; hyper1
%beta = 0.05; hyper2

GRID = {logspace(-4,-0.3,20),logspace(-4,-0.3,20)};

IND = find(triu(ones(CELLS,CELLS),1));

figure('position',[100         892        2102         487]);

subplot(1,5,1);imagesc(nodiag(real_corrmat));axis image;colorbar;
title('real')

[hypers_sample, bestDelta, visited, losses, quality_sample] = cove.crossEstimateHyper(DATA,size(DATA,1),'sample');%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat_sample,M_sample,V_sample,extras_sample,covmat_sample] = cove.estimate(DATA,size(DATA,1),'sample',[]); % alpha = 0.05; beta = 0.25
subplot(1,5,2);imagesc(nodiag(cormat_sample));axis image;colorbar;
title(sprintf('sample (%f)',corr(real_corrmat(IND),cormat_sample(IND))));

[hypers_diag, bestDelta_diag, visited_diag, losses_diag, quality_diag] = cove.crossEstimateHyper(DATA,size(DATA,1),'diag',GRID);%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat_diag,M_diag, V_diag, extras_diag,covmat_diag] = cove.estimate(DATA,size(DATA,1),'diag',hypers_diag); % alpha = 0.05; beta = 0.25
subplot(1,5,3);imagesc(nodiag(cormat_diag));axis image;colorbar;
title(sprintf('yatsenko diag (%f)',corr(real_corrmat(IND),cormat_diag(IND))));

[hypers_glasso, bestDelta_glasso, visited_glasso, losses_glasso, quality_glasso] = cove.crossEstimateHyper(DATA,size(DATA,1),'glasso',GRID(1));%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat_glasso,M_glasso,V_glasso,extras_glasso,covmat_glasso] = cove.estimate(DATA,size(DATA,1),'glasso',hypers_glasso); % alpha = 0.05; beta = 0.25
subplot(1,5,4);imagesc(nodiag(cormat_glasso));axis image;colorbar;
title(sprintf('yatsenko glasso (%f)',corr(real_corrmat(IND),cormat_glasso(IND))));

[hypers_lvglasso, bestDelta_lvglasso, visited_lvglasso, losses_lvglasso, quality_lvglasso] = cove.crossEstimateHyper(DATA,size(DATA,1),'lv-glasso',GRID);%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
[cormat_lvglasso,M_lvglasso,V_lvglasso,extras_lvglasso,covmat_lvglasso] = cove.estimate(DATA,size(DATA,1),'lv-glasso',hypers_lvglasso); % alpha = 0.05; beta = 0.25
subplot(1,5,5);imagesc(nodiag(cormat_lvglasso));axis image;colorbar;
title(sprintf('yatsenko model (%f)',corr(real_corrmat(IND),cormat_lvglasso(IND))));

figure('position',[100        100        2102         487]);
subplot(1,5,1);imagesc(nodiag(real_corrmat));axis image;colorbar;
title('real')

a = -corrcov(inv(covmat_diag));
subplot(1,5,2);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('sample sparse (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));

a = -corrcov(inv(covmat_diag));
subplot(1,5,3);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('yatsenko diag (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));

a = -corrcov(extras_glasso.S);
subplot(1,5,4);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('yatsenko glasso (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));

a = -corrcov(extras_lvglasso.S);
subplot(1,5,5);imagesc(nodiag(a));axis image;colorbar;
title(sprintf('yatsenko lv-glasso (%f, %f)',corr(real_corrmat(IND),a(IND)),edgeoverlap(real_corrmat,a,1)));



