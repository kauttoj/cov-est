
function [cormat,CELLS,label_str,testvals] = klab_compute_corrmat(dataOut,class_labels,THRESHOLD,BINNING,TYPE,CROSSVAL_ITERATIONS)

% Berkson's paradox

warning('on','all');

if nargin<6
    CROSSVAL_ITERATIONS = 0;
end

rng(666);

if nargin==0
    
    figure('position',[15         870        2322         628]);    
    
    [DATA,cormat_true,covmat_true] = make_random_network(100,5000,5); % nCells,nTrials,Bins
    %subplot(1,4,1);imagesc(nodiag(-corrcov(inv(covmat_true))));axis image;colorbar;
        
        
    [cormat3,M3, V3, extras3,covmat] = cove.estimate(DATA,size(DATA,1),'sample',[]); % alpha = 0.05; beta = 0.25
    covmat = (covmat + covmat')/2;
    cormat_partial3 = -corrcov(inv(covmat));
    subplot(1,3,1);imagesc(nodiag(cormat_partial3));axis image;colorbar;
    
    [hypers, bestDelta, visited, losses] = cove.crossEstimateHyper(DATA,size(DATA,1),'lv-glasso',{exp(-6:.1:-0.5),exp(-6:.1:-0.5)});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
    [cormat1,M1, V1, extras1,covmat] = cove.estimate(DATA,size(DATA,1),'lv-glasso',hypers); % alpha = 0.05; beta = 0.25
    covmat = (covmat + covmat')/2;
    cormat_partial1 = -corrcov(inv(covmat));
    subplot(1,3,2);imagesc(nodiag(cormat_partial1));axis image;colorbar;
    
    [hypers, bestDelta, visited, losses] = cove.crossEstimateHyper(DATA,size(DATA,1),'diag',{exp(-6:0.1:-0.5),exp(-6:0.1:-0.5)});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
    [cormat2,M2, V2, extras2,covmat] = cove.estimate(DATA,size(DATA,1),'diag',hypers); % alpha = 0.05; beta = 0.25
    covmat = (covmat + covmat')/2;
    cormat_partial2 = -corrcov(inv(covmat));
    subplot(1,3,3);imagesc(nodiag(cormat_partial2));axis image;colorbar;

    return;
    
end

DT = 1/dataOut.AquisitionParams.framerate;

N_bin_frames = round(BINNING/DT);

N_labels = length(dataOut.blockLabels);

frames = [];
for trial = 1:length(dataOut.trial)
   for repeat = 1:dataOut.trial(trial).stimuliRepeats       
       for stim = 1:N_labels
            frames(end+1) = sum((dataOut.trial(trial).stimulus_vector==stim).*(dataOut.trial(trial).repetition_vector==repeat));
       end
   end
end
N_blocks = median(frames);

bin_indices = 1:N_bin_frames:N_blocks;
bin_indices(end+1) = round(N_blocks);
bin_indices = unique(bin_indices);

N_bin_indices = diff(bin_indices);

%bin_window = mean(N_bin_indices)*DT;

if N_bin_indices(end)<mean(N_bin_indices(1:end-1)*0.5)
    N_bin_indices(end) = [];
    bin_indices(end) = [];
    warning('Removed last bin with less than 50% of frames!!')
end

N_bins = length(bin_indices)-1;

N_max_bad_frames = floor(N_bin_indices*0.25);

label = nan(1,N_labels);

if not(size(class_labels,1)==N_labels)
    error('Number of stimulus blocks does not match!')
end

k=0;
kk=0;
for i=1:N_labels
    for j=1:size(class_labels,1)
        if strcmp(class_labels{j,1},dataOut.blockLabels{i})
            k=k+1;
            if class_labels{j,2}>0
                kk=kk+1;
                label(i) = kk;
                label_str{i} = dataOut.blockLabels{i};
            else
                label(i) = 0;
            end
        end
    end
end
if k<N_labels
    error('not all labels were present in dataOut!');
end
valid = label>0;
label_ind = find(valid);
label_str = dataOut.blockLabels(label_ind);

N_labels = length(label_ind);

if ~isempty(THRESHOLD)
    if THRESHOLD==1
        CELLS = dataOut.stats.global.responsive_cells_p001_fdr_average(:)';
    else
        CELLS = THRESHOLD(:)';
    end
else
    CELLS = 1:dataOut.totalNumCells;
end
N_cells = length(CELLS);

DATA = zeros(N_cells,N_labels,N_bins,dataOut.totalNumTrials);
DATA_isbad = DATA;

k=0;
for trial = 1:length(dataOut.trial)
   for repeat = 1:dataOut.trial(trial).stimuliRepeats
       k=k+1;
       k1=0;
       for stim = label_ind
            k1=k1+1;
            vec = find((dataOut.trial(trial).stimulus_vector==stim).*(dataOut.trial(trial).repetition_vector==repeat));
            signal = dataOut.trial(trial).signal_deconv(CELLS,vec);            
            isbad = ~dataOut.trial(trial).isValidFrame(vec);
            for bin = 1:N_bins   
                binind = bin_indices(bin):(bin_indices(bin+1)-1);
                DATA(:,k1,bin,k) = mean(signal(:,binind),2);
                if sum(isbad(binind))>N_max_bad_frames(bin)                  
                    DATA_isbad(:,k1,bin,k) = 1;
                end
            end
       end
   end
end

grid1 = exp(-6:.05:-0.5);
grid2 = exp(-6:.05:-0.5);
    
ss = size(DATA);

if strcmp(TYPE,'sample')
    
    covmat = nan(N_cells,N_cells,N_labels);
    
    DATA_mean = mean(DATA,4);
    
    DATA = bsxfun(@minus,DATA,DATA_mean);
    
    assert(max(abs(flatten(mean(DATA,4))))<1e-13);
    
    DATA = reshape(DATA,[ss(1),ss(2),size(DATA,3)*size(DATA,4)]);
    
    for k=1:N_labels
        
        %ind = find(not(isnan(squeeze(DATA(1,k,:)))));
        
        data = squeeze(DATA(:,k,:))';
        
        assert(nnz(isnan(data))==0);
        
        fprintf('Label %i/%i with %i valid repeats (%i bad)\n',k,length(label_ind),size(data,1),ss(3)*ss(4)-size(data,1));
        
        covmat(:,:,k) = cov(data);
        
    end
    
    covmat = mean(covmat,3);
    
    cormat = corrcov(covmat);        
    
else
    
    %X - data: nBins * nConds * nTrials * nCells
    addpath('klab_noisecorrelation_tools');
    
    DATA1 = permute(DATA,[3,2,4,1]);    
    
    if strcmp(TYPE,'yatsenko_model')
        [hypers, bestDelta, visited, losses] = cove.crossEstimateHyper(DATA1,size(DATA1,1),'lv-glasso',{grid1,grid2});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
        [cormat,M, V, extras,covmat] = cove.estimate(DATA1,size(DATA1,1),'lv-glasso',hypers); % alpha = 0.05; beta = 0.25
    elseif strcmp(TYPE,'yatsenko_diag')
        [hypers, bestDelta, visited, losses] = cove.crossEstimateHyper(DATA1,size(DATA1,1),'diag',{grid1,grid2});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
        [cormat,M, V, extras,covmat] = cove.estimate(DATA1,size(DATA1,1),'diag',hypers); % alpha = 0.05; beta = 0.25
    elseif strcmp(TYPE,'yatsenko_sample')
        [cormat,M, V, extras,covmat] = cove.estimate(DATA1,size(DATA1,1),'sample',[]); % alpha = 0.05; beta = 0.25
    elseif strcmp(TYPE,'yatsenko_L1')        
        [hypers, bestDelta, visited, losses] = cove.crossEstimateHyper(DATA1,size(DATA1,1),'L1',{grid1});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
        [cormat,M, V, extras,covmat] = cove.estimate(DATA1,size(DATA1,1),'L1',hypers); % alpha = 0.05; beta = 0.25            
    else
        error('unknown method!');
    end
    
end

IND = find(triu(ones(N_cells,N_cells),1));

testvals=[];
if CROSSVAL_ITERATIONS>1
    testvals = do_crossval_iterations(CROSSVAL_ITERATIONS,DATA,N_cells,N_labels,TYPE,grid1,grid2,IND);
end

m = cormat(IND);
fprintf('Summary: %i cells, %i stimuli, %i trials (%i bins x %i repeats), mean correlation %f, 99%% [%f,%f]\n\n', ...
    ss(1),ss(2),ss(3)*ss(4),ss(3),ss(4),mean(m(:)),prctile(m,0.5),prctile(m,99.5));

%%
covmat_old = covmat;
covmat = (covmat + covmat')/2;

assert(max(abs(flatten(covmat-covmat_old)))<1e-10);

try
    cormat_partial = -corrcov(inv(covmat));
catch err
    warning('\nFailed to create partial correlation matrix!: %s\n',err.message);
    cormat_partial = nan;
end

if strcmp(TYPE,'yatsenko_model')
    cormat = {cormat,cormat_partial,-corrcov(extras.S)};
else
    cormat = {cormat,cormat_partial};
end

end

function mat = flatten(mat)
    mat = mat(:);
end

function mat = nodiag(mat)
    mat = mat - diag(diag(mat));
end

function testvals = do_crossval_iterations(CROSSVAL_ITERATIONS,DATA,N_cells,N_labels,TYPE,grid1,grid2,IND)

DATA_orig = DATA;

nTrials = size(DATA,4);
cv = cvpartition(nTrials,'KFold',2);

addpath('klab_noisecorrelation_tools');

for ITER = 1:CROSSVAL_ITERATIONS
    
    fprintf('-------- CV loop %i of %i\n',ITER,CROSSVAL_ITERATIONS);
    
    cv = cv.repartition;
    
    for part = 1:2
        
        if part == 1
            DATA = DATA_orig(:,:,:,cv.training(1));
        else
            DATA = DATA_orig(:,:,:,cv.test(1));
        end
        ss = size(DATA);
        
        if strcmp(TYPE,'sample')
            
            covmat = nan(N_cells,N_cells,N_labels);
            
            DATA_mean = mean(DATA,4);
            
            DATA = bsxfun(@minus,DATA,DATA_mean);
            
            assert(max(abs(flatten(mean(DATA,4))))<1e-13);
            
            DATA = reshape(DATA,[ss(1),ss(2),size(DATA,3)*size(DATA,4)]);
            
            for k=1:N_labels
                
                data = squeeze(DATA(:,k,:))';
                
                assert(nnz(isnan(data))==0);
                
                covmat(:,:,k) = cov(data);
                
            end
            
            covmat = mean(covmat,3);            
            cormat{1} = corrcov(covmat);
            cormat{2} = -corrcov(inv(covmat));
            
        else
            
            %X - data: nBins * nConds * nTrials * nCells            
            
            DATA1 = permute(DATA,[3,2,4,1]);
            
            if strcmp(TYPE,'yatsenko_model')
                hypers = cove.crossEstimateHyper(DATA1,size(DATA1,1),'lv-glasso',{grid1,grid2});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
                [cormat{part},~,~,extras{part},covmat{part}] = cove.estimate(DATA1,size(DATA1,1),'lv-glasso',hypers); % alpha = 0.05; beta = 0.25
            elseif strcmp(TYPE,'yatsenko_diag')
                hypers = cove.crossEstimateHyper(DATA1,size(DATA1,1),'diag',{grid1,grid2});%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
                [cormat{part},~,~,extras{part},covmat{part}] = cove.estimate(DATA1,size(DATA1,1),'diag',hypers); % alpha = 0.05; beta = 0.25
            elseif strcmp(TYPE,'yatsenko_sample')
                [cormat{part},~,~,extras{part},covmat{part}] = cove.estimate(DATA1,size(DATA1,1),'sample',[]); % alpha = 0.05; beta = 0.25
            else
               error('unknown type!') 
            end
            
        end
        
    end
    testvals.correlation.pearson(ITER) = corr(cormat{1}(IND),cormat{2}(IND),'type','pearson');
    testvals.correlation.spearman(ITER) = corr(cormat{1}(IND),cormat{2}(IND),'type','spearman');
    
    net1 = double(abs(cormat{1})>prctile(abs(cormat{1}(IND)),95));
    net2 = double(abs(cormat{2})>prctile(abs(cormat{2}(IND)),95));
    
    testvals.correlation.top5percent(ITER) = sum(sum(net1(IND).*net2(IND)))/(0.5*sum(sum(net1(IND)))+0.5*sum(sum(net1(IND))));
    
    a = -corrcov(inv(covmat{1}));
    b = -corrcov(inv(covmat{2}));
    testvals.partial_correlation.pearson(ITER) = corr(a(IND),b(IND),'type','pearson');
    testvals.partial_correlation.spearman(ITER) = corr(a(IND),b(IND),'type','spearman');
    
    net1 = double(abs(a)>prctile(abs(a(IND)),95));
    net2 = double(abs(b)>prctile(abs(b(IND)),95));
    
    testvals.partial_correlation.top5percent(ITER) = sum(sum(net1.*net2))/(0.5*sum(sum(net1))+0.5*sum(sum(net1)));
    
end

end