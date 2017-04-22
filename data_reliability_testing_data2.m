% simulation testbed

clearvars
close all;
clc;

FILENAME = '2093_NC_170414';

%A = load('data\traces-pre5-0001.mat');
%DATA = A.tuples(23).trace_segments;
%DATA = DATA(1:8,[1,4,5],:,:);

A=load('data\2093_NC_170414_data.mat');
%A=load('data\2093_NC_170408_data.mat');
DATA = permute(A.DATA,[3,2,4,1]);

DATA_orig = DATA + (1e-6)*rand(size(DATA)); % just in case there are null responses in one fold

test_params{1} = {10,200,1};
test_params{2} = {10,200,2};
test_params{3} = {10,200,3:4};
test_params{4} = {10,200,1:4};

test_methods = {'sample','diag','lv-glasso'};

test_results = [];

%X - data: nBins * nConds * nTrials * nCells

%alpha = 0.002; hyper1
%beta = 0.05; hyper2        

GRID = {logspace(-3.2,-0.1,20),logspace(-3.2,-0.1,20)};

all_seeds = randi(100000,1,100);

for TEST_LOOP = 1:length(test_params)    
    
    CELLS = test_params{TEST_LOOP}{2};
    CONDITIONS = test_params{TEST_LOOP}{3};
    
    for ITER = 1:test_params{TEST_LOOP}{1}
        
        rng(all_seeds(ITER));
        
        try
            
            cells_sample = randsample(size(DATA_orig,4),min(size(DATA_orig,4),CELLS));
            
            test_results{TEST_LOOP,ITER}.cells_sample = cells_sample;
            
            DATA = DATA_orig(:,CONDITIONS,:,cells_sample);
            
            %%%
            
            for k = 1:length(test_methods)
                
                [hypers, bestDelta, visited, losses, qualitytest] = cove.crossEstimateHyper(DATA,size(DATA,1),test_methods{k},GRID);%{10.^(-linspace(0.2,5,15)),10.^(-linspace(0.05,5,25))});
                [cormat,~,~,extras,covmat] = cove.estimate(DATA,size(DATA,1),test_methods{k},hypers); % alpha = 0.05; beta = 0.25
                
                test_results{TEST_LOOP,ITER}.data{k}.method = test_methods{k};
                test_results{TEST_LOOP,ITER}.data{k}.hypers = hypers;
                test_results{TEST_LOOP,ITER}.data{k}.qualitytest = qualitytest;
                test_results{TEST_LOOP,ITER}.data{k}.cormat = cormat;
                test_results{TEST_LOOP,ITER}.data{k}.extras = extras;
                test_results{TEST_LOOP,ITER}.data{k}.covmat = covmat;
                
            end
            
        catch err
            warning('failed loop!');
        end
        
    end
    
    save([FILENAME,'_test_results.mat'],'test_results','-v7.3');
    
end



