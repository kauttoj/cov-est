
% A=load('data\2093_NC_170408_data.mat');
% DATA = permute(A.DATA,[3,2,4,1]);

FILENAME = 'yatsenko_data';
A = load('data\traces-pre5-0001.mat');
DATA = A.tuples(18).trace_segments;
DATA = double(DATA(1:8,[1,2,4,5],:,:));

DATA_orig = DATA + (1e-10)*rand(size(DATA));

CONDITIONS = 1:4;
ITERATIONS = 50;

cond_estimate = zeros(30,30);

CELLS_arr = round(linspace(20,250,30));
TRIALS_arr = fliplr(round(linspace(10,143,30)));

k1=0;
for CELLS = CELLS_arr
    k1=k1+1;
    
    k2=0;
    for TRIALS = TRIALS_arr
        k2=k2+1;
    
        c = nan(1,ITERATIONS);
        
        for iter = 1:ITERATIONS
        
                cells_sample = randsample(size(DATA_orig,4),min(size(DATA_orig,4),CELLS));    
                trials_sample = randsample(size(DATA_orig,3),min(size(DATA_orig,4),TRIALS));
            
                DATA = DATA_orig(:,CONDITIONS,trials_sample,cells_sample);

                [cormat,~,~,extras,covmat] = cove.estimate(DATA,size(DATA,1),'sample',[]); % alpha = 0.05; beta = 0.25
                                
                c(iter) = extras.Cond;
                
        end
        
        cond_estimate(k1,k2) = median(c);
        
    end
end

a = cond_estimate/min(min(cond_estimate));
a(a>1000)=nan;

figure('position',[ 425         829        1246         669]);
pcolor(log10(a));
set(gca,'ytick',1:length(CELLS_arr),'yticklabel',CELLS_arr,'xtick',1:length(TRIALS_arr),'xticklabel',TRIALS_arr*8);
xlabel 'Trials (with 6 bins)'
ylabel 'Neurons'
colormap(jet(500))
b = colorbar();
c = get(b,'ytick');
set(b,'yticklabel',round(10.^c));
title('Yatsenko (all 4 conditions), normalized cov. matrix condition number (median of 50 iterations)','interpreter','none')




