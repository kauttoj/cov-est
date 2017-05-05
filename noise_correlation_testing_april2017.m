clc;

clearvars;
close all;

%addpath('C:\Users\skOpti_6\Documents\MATLAB\klab_cov-est'); no 

BIN_WIDTH = 0.250; % in seconds
CORRECT_DATA = 2; % 0 = no correction, 1 = use neuropil, 2 = use signal baseline
SHUFFLE_TRIALS = 1; % if 1, file and repeat identities are shuffled

% DATASETS
ROOT_PATH1 = 'Z:\ProcessedDataArchive\noise_correlations\2093_NC_170420_suite2p_processed\processed_suite2p\analysis';
DATAOUT_file1 = [ROOT_PATH1,filesep,'2093_NC_170420_PINKNOISE_dataOut.mat'];
A1 = load(DATAOUT_file1);
[~,dataname1] = fileparts(DATAOUT_file1);

ROOT_PATH2 = 'Z:\ProcessedDataArchive\noise_correlations\2093_NC_170421_suite2p_processed\processed_suite2p\analysis';
DATAOUT_file2 = [ROOT_PATH2,filesep,'2093_NC_170421_PINKNOISE_dataOut.mat'];
A2 = load(DATAOUT_file2);
[~,dataname2] = fileparts(DATAOUT_file2);

datanames = {dataname1(1:end-18),dataname2(1:end-18)};

% STIMULI
class_labels = {...
    'forest_biking2',...
    'manhattan1',...
    'pink_noise_rate0.80_vid1',...
    'pink_noise_rate0.80_vid2',...
    };

%%

% PARSE DATA
DATA1 = klab_get_data(A1.dataOut,class_labels,1,BIN_WIDTH);

DATA2 = klab_get_data(A2.dataOut,class_labels,1,BIN_WIDTH);

if SHUFFLE_TRIALS    
    [rep,fil] = meshgrid(1:size(DATA1.DATA,4),1:size(DATA1.DATA,2));
    ind = randperm(numel(fil));
    fil(:) = fil(ind);
    rep(:) = rep(ind);
    
    DATA1_old = DATA1;
    
    for i=1:size(fil,1)
        for j=1:size(rep,2)
            DATA1.DATA(:,i,:,j,:) = DATA1_old.DATA(:,fil(i,j),:,rep(i,j),:);
            DATA1.DATA_neuropil(:,i,:,j,:) = DATA1_old.DATA_neuropil(:,fil(i,j),:,rep(i,j),:);
            DATA1.DATA_baseline(:,i,:,j,:) = DATA1_old.DATA_baseline(:,fil(i,j),:,rep(i,j),:);
            DATA1.DATA_isbad(:,i,:,j,:) = DATA1_old.DATA_isbad(:,fil(i,j),:,rep(i,j),:);            
        end
    end
    
    assert(max(abs(sort(DATA1_old.DATA(:))-sort(DATA1.DATA(:))))<eps);
    
    [rep,fil] = meshgrid(1:size(DATA2.DATA,4),1:size(DATA2.DATA,2));
    ind = randperm(numel(fil));
    fil(:) = fil(ind);
    rep(:) = rep(ind);
    
    DATA2_old = DATA2;
    
    for i=1:size(fil,1)
        for j=1:size(rep,2)           
            DATA2.DATA(:,i,:,j,:) = DATA2_old.DATA(:,fil(i,j),:,rep(i,j),:);
            DATA2.DATA_neuropil(:,i,:,j,:) = DATA2_old.DATA_neuropil(:,fil(i,j),:,rep(i,j),:);
            DATA2.DATA_baseline(:,i,:,j,:) = DATA2_old.DATA_baseline(:,fil(i,j),:,rep(i,j),:);
            DATA2.DATA_isbad(:,i,:,j,:) = DATA2_old.DATA_isbad(:,fil(i,j),:,rep(i,j),:);
        end
    end    
    
    assert(max(abs(sort(DATA2_old.DATA(:))-sort(DATA2.DATA(:))))<eps);
    
end


%% TRY TO RESCALE DATA TO COMPENSATE WATER-LOSS

if CORRECT_DATA>0
    
    if CORRECT_DATA==1
        correction_data1 = DATA1.DATA_neuropil;
        correction_data2 = DATA2.DATA_neuropil;
    elseif CORRECT_DATA==2
        correction_data1 = DATA1.DATA_baseline;
        correction_data2 = DATA2.DATA_baseline;
    else
        error('unknown');
    end
    
    allmeans1 = [];
    for i=1:size(DATA1.DATA,1)
        a = squeeze(correction_data1(i,:,:,:,:));
        m = mean(a(:));
        oldvar = var(vectorize(DATA1.DATA(i,:,:,:,:)));
        coef = (m./correction_data1(i,:,:,:,:));
        DATA1.DATA(i,:,:,:,:) = DATA1.DATA(i,:,:,:,:).*coef;
        assert(nnz(isnan(DATA1.DATA(i,:,:,:,:)))==0);
        newvar = var(vectorize(DATA1.DATA(i,:,:,:,:)));
        allmeans1(i,1:5) = [min(coef(:)),mean(coef(:)),max(coef(:)),oldvar,newvar];
    end
    
    allmeans2 = [];
    for i=1:size(DATA2.DATA,1)
        a = squeeze(correction_data2(i,:,:,:,:));
        m = mean(a(:));
        oldvar = var(vectorize(DATA2.DATA(i,:,:,:,:)));
        coef = (m./correction_data2(i,:,:,:,:));
        DATA2.DATA(i,:,:,:,:) = DATA2.DATA(i,:,:,:,:).*coef;
        assert(nnz(isnan(DATA2.DATA(i,:,:,:,:)))==0);
        newvar = var(vectorize(DATA2.DATA(i,:,:,:,:)));
        allmeans2(i,1:5) = [min(coef(:)),mean(coef(:)),max(coef(:)),oldvar,newvar];
    end
end

DATA1_orig = DATA1;
DATA2_orig = DATA2;

N_cells_orig = [size(DATA1.DATA,1),size(DATA2.DATA,1)];

%% GET SHARED CELLS
C = load('Z:\ProcessedDataArchive\noise_correlations\2093_NC_170421_suite2p_processed\processed_suite2p\registration 20170424T114647\2093_NC_170421_globalCellID.mat');

mat = C.globalID;

%% remove double matches
mat((mat(:,2)==67),:)=[];

assert(length(unique(mat(:,1))) + length(unique(mat(:,1))) == 2*size(mat,1));
%% only keep data of common cells

mat(~ismember(mat(:,1),DATA1.CELL_IDs),1)=0;
mat(~ismember(mat(:,2),DATA2.CELL_IDs),2)=0;

mat((sum(mat>0,2)<2),:)=[];

ind1=[];
ind2=[];
for i=1:size(mat,1)
    ind1(end+1) = find(DATA1.CELL_IDs==mat(i,1));
    ind2(end+1) = find(DATA2.CELL_IDs==mat(i,2));
end
DATA1.DATA = DATA1.DATA(ind1,:,:,:,:);
DATA2.DATA = DATA2.DATA(ind2,:,:,:,:);
DATA1.CELL_IDs = DATA1.CELL_IDs(ind1);
DATA2.CELL_IDs = DATA2.CELL_IDs(ind2);
DATA1.DATA_isbad = DATA1.DATA_isbad(ind1,:,:,:,:);
DATA2.DATA_isbad = DATA2.DATA_isbad(ind2,:,:,:,:);

N_cells = size(DATA1.DATA,1);

%% START TESTS
rng(666); % fixed seed

stimuli = {1,2,3,4,1:4}; % choose stimuli
ITERATIONS = 10; % cv runs
NULL_ITERATIONS = 150; % randomized neurons ids
cv{1} = cvpartition(size(DATA1.DATA,2),'KFold',2);
cv{2} = cvpartition(size(DATA2.DATA,2),'KFold',2);
density_sweep = 50:99; % percentage for binarization
results_cross_session = struct();
results_inside_session = cell(2,ITERATIONS);

IND_common = find(triu(ones(N_cells),1));
IND{1} = find(triu(ones(N_cells_orig(1)),1));
IND{2} = find(triu(ones(N_cells_orig(2)),1));

% cross session
k=0;
for stim = stimuli
    k=k+1;
    
    [cormat1,part_cormat1] = compute_correlations(DATA1.DATA(:,:,stim{1},:,:));
    [cormat2,part_cormat2] = compute_correlations(DATA2.DATA(:,:,stim{1},:,:));
    
    vals1 = squareform(nodiag(cormat1));
    vals2 = squareform(nodiag(cormat2));
            
    part_vals1 = squareform(nodiag(part_cormat1));
    part_vals2 = squareform(nodiag(part_cormat2));    
    
    density_corr = arrayfun(@(x) corr2(vals1.*(vals1>prctile(vals1,x)),vals2.*(vals2>prctile(vals2,x))),density_sweep);
    density_part_corr = arrayfun(@(x) corr2(part_vals1.*(part_vals1>prctile(part_vals1,x)),vals2.*(part_vals2>prctile(part_vals2,x))),density_sweep);   
    
    results_cross_session(k).stim = stim{1};
    
    results_cross_session(k).corrvals1 = vals1;
    results_cross_session(k).corrvals2 = vals2;    
    results_cross_session(k).similarity_corr = corr2(vals1,vals2);
    results_cross_session(k).similarity_partial_corr = corr2(part_vals1,part_vals2);
    results_cross_session(k).density_corr = density_corr;
    results_cross_session(k).density_part_corr = density_part_corr;    
    results_cross_session(k).differences_corr = vals1-vals2;
    results_cross_session(k).differences_part_corr  = part_vals1-part_vals2;
    
    for null_iter = 1:NULL_ITERATIONS
        
        perm = randperm(N_cells);
        
        cormat1_rand = cormat1(perm,perm);
        
        part_cormat1_rand = part_cormat1(perm,perm);
        
        vals1 = cormat1_rand(IND_common);
        vals2 = cormat2(IND_common);
                
        part_vals1 = part_cormat1_rand(IND_common);
        part_vals2 = part_cormat2(IND_common);
        
        density_corr = arrayfun(@(x) corr2(vals1.*(vals1>prctile(vals1,x)),vals2.*(vals2>prctile(vals2,x))),density_sweep);
        density_part_corr = arrayfun(@(x) corr2(part_vals1.*(part_vals1>prctile(part_vals1,x)),vals2.*(part_vals2>prctile(part_vals2,x))),density_sweep);        
                
        results_cross_session(k).similarity_corr_nullvals(null_iter) = corr2(vals1,vals2);
        results_cross_session(k).similarity_partial_corr_nullvals(null_iter) = corr2(part_vals1,part_vals2);
        results_cross_session(k).density_corr_nullvals{null_iter} = density_corr;
        results_cross_session(k).density_part_corr_nullvals{null_iter} = density_part_corr;           
        results_cross_session(k).differences_corr_nullvals{null_iter} = vals1-vals2;
        results_cross_session(k).differences_part_corr_nullvals{null_iter}  = part_vals1-part_vals2;
        
    end
    
end

% inside session
for dataset = 1:2
    
    if dataset==1
        DATA = DATA1_orig.DATA;
    else
        DATA = DATA2_orig.DATA;
    end
    
    k=0;
    for stim = stimuli

        k=k+1;
        
        fprintf('starting: dataset %i, stim %i\n',dataset,k);
                
        for cv_iter = 1:ITERATIONS
            
%             if cv_iter==1
%                 files1 = 1:sum(cv{dataset}.training(1));
%                 files2 = length(files1) + (1:sum(cv{dataset}.training(2)));
%             else
                files1 = cv{dataset}.training(1);
                files2 = cv{dataset}.training(2);
            %end
            
            [cormat1,part_cormat1] = compute_correlations(DATA(:,files1,stim{1},:,:));
            [cormat2,part_cormat2] = compute_correlations(DATA(:,files2,stim{1},:,:));
            
            vals1 = cormat1(IND{dataset});
            vals2 = cormat2(IND{dataset});
            
            part_vals1 = part_cormat1(IND{dataset});
            part_vals2 = part_cormat2(IND{dataset});
                        
            results_inside_session{dataset,cv_iter}(k).stim = stim{1};
            
            results_inside_session{dataset,cv_iter}(k).corrvals1 = vals1;
            results_inside_session{dataset,cv_iter}(k).corrvals2 = vals2;                        
            
            results_inside_session{dataset,cv_iter}(k).similarity_corr = corr2(vals1,vals2);
            results_inside_session{dataset,cv_iter}(k).similarity_partial_corr = corr2(part_vals1,part_vals2);
            results_inside_session{dataset,cv_iter}(k).differences_corr = vals1-vals2;
            results_inside_session{dataset,cv_iter}(k).differences_part_corr  = part_vals1-part_vals2;
           
            cv{dataset} = cv{dataset}.repartition;
            
            for null_iter = 1:NULL_ITERATIONS
                
                perm = randperm(size(DATA,1));                
                
                cormat1_rand = cormat1(perm,perm);
                
                part_cormat1_rand = part_cormat1(perm,perm);
                
                vals1 = cormat1_rand(IND{dataset});
                vals2 = cormat2(IND{dataset});
                
                part_vals1 = part_cormat1_rand(IND{dataset});
                part_vals2 = part_cormat2(IND{dataset});
                
                results_inside_session{dataset,cv_iter}(k).similarity_corr_nullvals(null_iter) = corr2(vals1,vals2);
                results_inside_session{dataset,cv_iter}(k).similarity_partial_corr_nullvals(null_iter) = corr2(part_vals1,part_vals2);
                results_inside_session{dataset,cv_iter}(k).differences_corr_nullvals{null_iter} = vals1-vals2;
                results_inside_session{dataset,cv_iter}(k).differences_part_corr_nullvals{null_iter}  = part_vals1-part_vals2;
                
            end            
            
        end
        
    end
end

%% plots
COLORS = [...
    0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];

mkdir('inside_session');
mkdir('between_sessions');

%% cross session

close all
fig1 = figure('position',[633   806   927   692]);
hold on;
title(sprintf('Correlation matrix similarity between sessions (%i neurons)',N_cells));


fig2 = figure('position',[633   806   927   692]);
hold on;
title(sprintf('Correlation differences (%i neurons)',N_cells));

leg1 = [];
leg2 = [];
leg_handle1 = [];
leg_handle2 = [];

for stim=1:length(stimuli)
    figure(fig1);
    differences_corr_nullvals = [];
    a = histogram(results_cross_session(stim).similarity_corr_nullvals,'normalization','probability','DisplayStyle','stairs','linewidth',2,'EdgeColor',COLORS(stim,:));
    axis tight;
    y=ylim();
    hold on;
    leg_handle1(end+1) = plot(results_cross_session(stim).similarity_corr*[1,1],y(2)*[0,1]/2,'linewidth',2,'Color',COLORS(stim,:));    
    leg1{end+1}=['stim: ',num2str(stimuli{stim})];

    figure(fig2);
    dist_real = abs(results_cross_session(stim).differences_corr)/sqrt(2);
    dist_nullvals = [];
    for null_iter=1:length(results_cross_session(stim).differences_corr_nullvals)
        dist_nullvals = [dist_nullvals;results_cross_session(stim).differences_corr_nullvals{null_iter}];
    end
    dist_nullvals = abs(dist_nullvals)/sqrt(2);  
    a = histogram(dist_nullvals,0:0.02:0.6,'normalization','probability','DisplayStyle','stairs','linewidth',2,'EdgeColor',COLORS(stim,:),'linestyle','--');
    leg_handle2(end+1) = histogram(dist_real,0:0.02:0.6,'normalization','probability','DisplayStyle','stairs','linewidth',3,'EdgeColor',COLORS(stim,:));           
    [h,p] = kstest2(dist_nullvals,dist_real);
    leg2{end+1}=sprintf('stim: %s (p=%1.2e)',num2str(stimuli{stim}),p);
end
figure(fig1);
x = xlim();
xlim([x(1),x(2)+range(x)*0.05]);
legend(leg_handle1,leg1,'location','best');
box on;
xlabel('Network similarity (upper triagonal)')
ylabel('Count density')
set(gca,'FontSize',14);
saveas(fig1,['between_sessions\corrmat_similarity.fig']);
saveas(fig1,['between_sessions\corrmat_similarity.png']);

figure(fig2);
legend(leg_handle2,leg2,'location','best');
box on;
set(gca,'yscale','log');
xlabel('Distance from x==y')

ylabel('Count density (log scale)')
set(gca,'FontSize',14);

saveas(fig2,['between_sessions\corrval_differences.fig']);
saveas(fig2,['between_sessions\corrval_differences.png']);

% scatterplots

k=0;
for stim=stimuli
    k=k+1;
    figure('position',[750   920   810   578]);
    scatter(results_cross_session(k).corrvals1,results_cross_session(k).corrvals2);
    axis tight;    
    x=xlim();y=ylim();    
    hold on;
    plot([-1,1],[-1,1],'k');
    xlim(x);ylim(y);
    title(sprintf('stim = %s, corr=%1.3f, %i neurons, %i edges',num2str(stim{1}),corr2(results_cross_session(k).corrvals1,results_cross_session(k).corrvals2),N_cells,length(results_cross_session(k).corrvals1))); 
    xlabel(datanames{1});
    ylabel(datanames{2});
    axis equal;
    box on;
    saveas(gcf,sprintf('between_sessions\\scatterplot_stimset%i.fig',k));
    saveas(gcf,sprintf('between_sessions\\scatterplot_stimset%i.png',k));
end

%% inside session

close all

for dataset = 1:2
    
    fig1 = figure('position',[633   806   927   692]);
    hold on;
    title(sprintf('%s, correlation matrix similarity between folds (%i neurons)',datanames{dataset},N_cells_orig(dataset)));
        
    fig2 = figure('position',[633   806   927   692]);
    hold on;
    title(sprintf('%s, correlation differences (%i neurons)',datanames{dataset},N_cells_orig(dataset)));
    
    leg1 = [];
    leg2 = [];
    leg_handle1 = [];
    leg_handle2 = [];
    
    for stim=1:length(stimuli)
        figure(fig1);
        similarity_corr = [];
        similarity_corr_nullvals=[];
        for i=1:size(results_inside_session,2)
            similarity_corr = [similarity_corr;results_inside_session{dataset,i}(stim).similarity_corr]; 
            similarity_corr_nullvals = [similarity_corr_nullvals,results_inside_session{dataset,i}(stim).similarity_corr_nullvals];
        end
        a = histogram(similarity_corr_nullvals,'normalization','probability','DisplayStyle','stairs','linewidth',2,'EdgeColor',COLORS(stim,:));
        axis tight;
        y=ylim();
        hold on;
        leg_handle1(end+1) = plot(mean(similarity_corr)*[1,1],y(2)*[0,1]/2,'linewidth',2,'Color',COLORS(stim,:));
        leg1{end+1}=['stim: ',num2str(stimuli{stim})];
        
        figure(fig2);
        dist_nullvals = [];
        dist_real = 0;
        for i=1:size(results_inside_session,2)
            dist_real = [dist_real;abs(results_inside_session{dataset,i}(stim).differences_corr)/sqrt(2)];
            for null_iter=1:10
                dist_nullvals = [dist_nullvals;abs(results_inside_session{dataset,i}(stim).differences_corr_nullvals{null_iter})/sqrt(2)];
            end
        end
        a = histogram(dist_nullvals,0:0.02:0.6,'normalization','probability','DisplayStyle','stairs','linewidth',2,'EdgeColor',COLORS(stim,:),'linestyle','--');
        leg_handle2(end+1) = histogram(dist_real,0:0.02:0.6,'normalization','probability','DisplayStyle','stairs','linewidth',3,'EdgeColor',COLORS(stim,:));
        [h,p] = kstest2(dist_nullvals,dist_real);
        leg2{end+1}=sprintf('stim: %s (p=%1.2e)',num2str(stimuli{stim}),p);
    end
    figure(fig1);
    x = xlim();
    xlim([x(1),x(2)+range(x)*0.05]);
    legend(leg_handle1,leg1,'location','best');
    box on;
    xlabel('Network similarity (upper triagonal)')
    ylabel('Count density')
    set(gca,'FontSize',14);
    saveas(fig1,['inside_session\',datanames{dataset},'_corrmat_similarity.fig']);
    saveas(fig1,['inside_session\',datanames{dataset},'_corrmat_similarity.png']);
    
    figure(fig2);
    legend(leg_handle2,leg2,'location','best');
    box on;
    set(gca,'yscale','log');
    xlabel('Distance from x==y')    
    ylabel('Count density (log scale)')
    set(gca,'FontSize',14);
    saveas(fig2,['inside_session\',datanames{dataset},'_corrval_differences.fig']);
    saveas(fig2,['inside_session\',datanames{dataset},'_corrval_differences.png']);    
end

% scatterplots

for dataset = 1:2
    for cv_iter = 1:2 % just first two
        k=0;
        for stim=stimuli
            k=k+1;
            figure('position',[750   920   810   578]);
            scatter(results_inside_session{dataset,cv_iter}(k).corrvals1,results_inside_session{dataset,cv_iter}(k).corrvals2);
            axis tight;
            x=xlim();y=ylim();
            hold on;
            plot([-1,1],[-1,1],'k');
            xlim(x);ylim(y);
            title(sprintf('stim = %s, corr=%1.3f, %i neurons, %i edges',num2str(stim{1}),corr2(results_inside_session{dataset,...
                cv_iter}(k).corrvals1,results_inside_session{dataset,cv_iter}(k).corrvals2),N_cells_orig(dataset),...
                length(results_inside_session{dataset,cv_iter}(k).corrvals1)));
            xlabel('Fold 1');
            ylabel('Fold 2');
            axis equal;
            box on;
            saveas(gcf,sprintf('inside_session/%s_scatterplot_stimset%i_cviter%i.fig',datanames{dataset},k,cv_iter));
            saveas(gcf,sprintf('inside_session/%s_scatterplot_stimset%i_cviter%i.png',datanames{dataset},k,cv_iter));
        end
    end
end
