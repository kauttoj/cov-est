function [hypers, bestDelta, visited, losses, extras] = crossEstimateHyper(X, evokedBins, reg, searchSpace)
% find the values of hyperparamters that minimize the cross-validated loss
% by K-fold cross-validation.
%
% INPUTS:
%     X = nBins * nDirs * nTrials * nCells
%     evoked bins - number of bins in each trial that are driven by stimulus
%     reg - structure specifying regularization of means, variances, and correlations
%     searchSpace - lists of valid values for each hyperparameter, in sequence
%
% OUTPUTS:
%     hypers - optimal values of hyperparameters
%     visited - the indices of visited hyperparameter values
%     losses  - the loss function values for these hyperparameter values

visited = [];  % visited indices
losses =  [];  % measured losses at visited indices
correlations = [];
edge_overlaps = [];
Rps = [];

K = 10;

if isempty(searchSpace) || strcmp(reg,'sample')
    [XTest,R,M,V] = arrayfun(@(k) estimate_([],k,K), 1:K,  'UniformOutput', false);
    [losses_arr,corrvals,overlaps,Rp] = cellfun(@(XTest,R,M,V) cove.vloss(XTest,R,M,V,0,reg,[]), XTest, R, M, V, 'UniformOutput', false);
    
    losses_arr = cell2mat(losses_arr);
    corrvals = cell2mat(corrvals);
    overlaps = cell2mat(overlaps);
    
    losses = mean(losses_arr);
    extras.correlation = corrvals;
    extras.edge_overlap = overlaps;
    
    extras.Rp = Rp;
    hypers = [];
    bestDelta = nan;
    return;
end

switch reg
    case 'diag'
        % 2 hypers:  variance shrinkage,  correlation shrinkage
        assert(length(searchSpace)==2, 'invalid hyperparameters')       
    case 'factor'
        % 2 hypers: variance shrink toward median, nlatent
        assert(length(searchSpace)==2, 'invalid hyperparameters') 
    case 'glasso'
        searchSpace = searchSpace(1);
        assert(length(searchSpace)==1, 'invalid hyperparameters') 
    case 'L1'        
        searchSpace = searchSpace(1);
        assert(length(searchSpace)==1, 'invalid hyperparameters') 
    case 'QUIC'
        searchSpace = searchSpace(1);
        assert(length(searchSpace)==1, 'invalid hyperparameters')         
    case 'lv-glasso'
        assert(length(searchSpace)==2, 'invalid hyperparameters')         
    otherwise
        error 'unknown covariance estimator'
end


dims = cellfun(@length, searchSpace);
nHypers = length(dims);
assert(nHypers>0)


[X1,Fval,Exitflag,Output] = patternsearch(objectiveFcn,X0,Aineq,Bineq,Aeq,Beq);

% % coarse random search to seed the local search
% fprintf 'random search: '
% decimate = 5;
% nRandom = ceil(prod(max(1,dims/decimate)/sum(dims>1)));
% fprintf('%d points\n', nRandom)
% ix = arrayfun(@(~) arrayfun(@(d) randi(d), dims), 1:nRandom, 'uni', false);
% cellfun(@visit,ix);
% 
% % pattern search for the optimum hyperparameter values
% disp 'pattern search'
% pattern = dims>1;
% step = min(floor(dims/2),ceil(dims/decimate*2));
% bestLoss = min(losses);
% while true
%     lastBestLoss = bestLoss;
%     [bestLoss,j] = min(losses);
%     assert(~isnan(bestLoss))
%     if bestLoss == lastBestLoss && all(step<=1)
%         break
%     end
%     step = ceil(step/2);
%     [best{1:nHypers}] = ind2sub(dims,visited(j));
%     ix = [best{:}];
%     % visit all nodes 1 step away in every direction
%     for b = 0:2^sum(pattern)-1
%         s = step;
%         s(pattern) = s(pattern).*(2*bitget(b,1:sum(pattern))-1);
%         visit(ix+s)
%     end
% end
% [~,j] = min(losses);
% [indices{1:nHypers}] = ind2sub(dims,visited(j));
% hypers = cellfun(@(h,ix) h(ix), searchSpace, indices);

extras.correlation = correlations{j};
extras.edge_overlap = edge_overlaps{j};
extras.Rp = Rps{j};

fprintf('final hyperparameters: %s\n', sprintf(' %g', hypers))

    function visit(ix)
        % visit a point on the search space specified by ix
        % but return if already visited
        
        ix = num2cell(max(1,min(dims,ix)));
        hypers_ = cellfun(@(x,i) x(i), searchSpace, ix);
        ix = sub2ind([dims ones(1,2-length(dims))],ix{:});
        alreadyVisited = ismember(ix,visited);
        if ~alreadyVisited
            fprintf('%2.4g  ', hypers_)
            [XTest,R,M,V] = arrayfun(@(k) estimate_(hypers_,k,K), 1:K,  'UniformOutput', false);
            visited(end+1) = ix;
            if size(X,2)==1 % && size(X,1)==1
                delta = 0;
            else
                % if multiple conditions, regularize variance estimate
                deltas = cellfun(@(XTest,R,M,V) cove.findBestDelta(XTest, R, M, V), XTest, R, M, V);
                delta = mean(deltas);
            end
            [losses_arr,corrvals,overlaps,Rp] = cellfun(@(XTest,R,M,V) cove.vloss(XTest,R,M,V,delta,reg,hypers_), XTest, R, M, V,'UniformOutput', false);            
            
            losses_arr = cell2mat(losses_arr);
            corrvals = cell2mat(corrvals);
            overlaps = cell2mat(overlaps);            
            
            losses(end+1) = mean(losses_arr);
            correlations{end+1} = corrvals;
            edge_overlaps{end+1} = overlaps;
            Rps{end+1} = Rp;
            
            if losses(end)==min(losses)
                bestDelta = delta;
            end
            fprintf(':  mean loss %g\n', losses(end))
        end
    end

    function [XTest,R,M,V] = estimate_(hypers,k,K)
        % compute cross-validation loss
        [XTrain,XTest] = cove.splitTrials(X,k,K);
        [R, M, V] = cove.estimate(XTrain, evokedBins, reg, hypers,0);
    end
end