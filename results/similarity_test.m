clc;
clearvars
close all

load 2093_NC_170408_test_results.mat
% load 2093_NC_170414_test_results.mat

for row = 1:size(test_results,1)
    for col = 1:size(test_results,2)
        fprintf('row %i of %i, column %i of %i\n',row,size(test_results,1),col,size(test_results,2));
        for k=1:3
            
            cormat = nodiag(test_results{row,col}.data{k}.cormat);            
            %cormat = -corrcov(test_results{row,col}.data{k}.extras.S);            
            
            [mat,ind] = concat_cell(test_results{row,col}.data{k}.qualitytest.Rp,1);
            
            cormat = cormat(ind);
            th = prctile(abs(cormat),99);
            cormat = (abs(cormat)>th);
            
            p = nan(size(mat,1),1);
            for kk=1:size(mat,1)
                p(kk) = signtest(mat(kk,:));
            end
            p005_edges(row,col,k) = nnz(p<0.05)/length(p);
            
            top1_overlap_prc(row,col,k) = sum(cormat.*(p<0.05))/sum(cormat);
            
        end
    end
end

%%

figure('position',[591         338         969        1160])

for k=1:3
    subplot(3,1,k);
    imagesc(top1_overlap_prc(:,:,k))
    colorbar
    title(sprintf('2093_NC_170414, 10-fold CV p<0.05 (signtest) of top 1%%, model=%s',test_results{1,1}.data{k}.method),'interpreter','none')
    set(gca,'ytick',[1:4],'yticklabel',{'movie 1','movie 2','noise 1&2','all four'},'xtick',1:10)
    ylabel('Stimulus')
    xlabel('Random 200 neuron group')
    set(gca,'fontsize',16);
end


%%

p005_edges_diff = [];
p001_edges_diff = [];
pval=[];
overlap=[];
for col = 1:size(test_results,2)
    fprintf('column %i of %i\n',col,size(test_results,2));
    for k=1:3
        
        [mat1,ind] = concat_cell(test_results{1,col}.data{k}.qualitytest.Rp,1);
        [mat2,ind] = concat_cell(test_results{2,col}.data{k}.qualitytest.Rp,1);
        
        cormat1 = nodiag(test_results{1,col}.data{k}.cormat);
        cormat2 = nodiag(test_results{2,col}.data{k}.cormat);
        
        cormat1 = cormat1(ind);
        th = prctile(abs(cormat1),99);
        cormat1 = (abs(cormat1)>th);
        
        cormat2 = cormat2(ind);
        th = prctile(abs(cormat2),99);
        cormat2 = (abs(cormat2)>th);
        
        cormat = (cormat1 + cormat2 > 0);
        
        p = nan(size(mat1,1),1);
        for kk=1:size(mat1,1) 
            p(kk) = signrank(mat1(kk,:),mat2(kk,:));
        end
        p005_edges_diff(col,k) = nnz(p<0.05);
        p001_edges_diff(col,k) = nnz(p<0.01);
        
        overlap(col,k) = sum((p<0.05).*cormat)/sum(cormat);
        
        pval(col,k) = lottery_pval(cormat,p<0.05);
                
    end
end

title(sprintf('2093_NC_170408, 10-fold CV p<0.01 (signrank), movie 1 vs. movie 2'),'interpreter','none')
