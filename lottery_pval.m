function [pval,thresholds] = lottery_pval(vec1,vec2,PRINT)

    if nargin==0
        
        resmat = nan(20,7);
        edgemat = nan(20,7);
        
        threshold_arr = linspace(0.005,0.30,40);%[0.5,1,5,10,15,20,25]/100;
        nodes_arr = 68;%round(linspace(40,200,20));
        
        k1=0;
        for nodes = nodes_arr
            edges = nodes*(nodes-1)/2;
            vec = zeros(1,edges);            
            k1=k1+1;
            k2=0;
            for threshold = threshold_arr
                i = round(edges*threshold);
                if i==0
                    continue;
                end
                fprintf('\n..nodes %i (%i edges), threshold %f (%i edges)',nodes,edges,threshold,i);
                k2=k2+1;        
                assert(i<length(vec));
                vec(1:i)=1;
                [~,thresholds] = lottery_pval(vec,vec,0);
                resmat(k1,k2) = thresholds(3);
                edgemat(k1,k2) = edges;
            end
        end
        fprintf('\ndone!\n')
        
        figure('position',[1000         837        1500         661]);
        
        subplot(1,2,1);
        imagesc((resmat));
        colorbar;
        title('number of matching edges at p=0.001');
        set(gca,'XTick',1:size(resmat,2),'xticklabel',threshold_arr,'ytick',1:size(resmat,1),'yticklabel',nodes_arr);
        xlabel('total ratio of edges')
        ylabel('neurons')

        subplot(1,2,2);
        imagesc(resmat./edgemat);
        colorbar;
        title('density of matching edges at p=0.001');
        set(gca,'XTick',1:size(resmat,2),'xticklabel',threshold_arr,'ytick',1:size(resmat,1),'yticklabel',nodes_arr);        
        xlabel('total ratio of edges')
        ylabel('neurons')

    end
    
    if nargin<3
        PRINT = 1;
    end

    M = 20000;
    
    vec1 = vec1~=0;
    vec2 = vec2~=0;
    
    if sum(vec1) ==0 || sum(vec2)==0
       error('empty vectors!') ;
    end
    
    N = length(vec1);
    NN = sum(vec1);
        
    n = sum(vec1.*vec2);
    
    null = zeros(1,M);
    for i=1:M
        vec1 = 0*vec1;
        vec1(randsample(N,NN,'false'))=1;
        null(i) = sum(vec1.*vec2);
    end
    
    pval = nnz(null>=n)/length(null);
    
    thresholds = ceil(prctile(null,[95,99,99.9]));
    
    if PRINT
        fprintf('pval = %f (p=0.01 threshold is %i)\n',pval,ceil(prctile(null,99)));
    end

end