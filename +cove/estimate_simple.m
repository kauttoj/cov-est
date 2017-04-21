function C = estimate_simple(C,reg,hypers)

switch reg
    case 'sample'
        assert(isempty(hypers),'invalid hyperparameters')
        % do nothing
        
    case 'diag'
        % 2 hypers:  variance shrinkage,  correlation shrinkage
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        r = (1-hypers(2))*corrcov(C) + hypers(2)*eye(size(C));
        C = v*r*v;
        
    case 'factor'
        % 2 hypers: variance shrink toward median, nlatent
        assert(length(hypers)==2, 'invalid hyperparameters')
        v = diag(C);
        v = diag(sqrt((1-hypers(1))*v + hypers(1)*median(v)));
        C = v*corrcov(C)*v;
        [L,psi] = cove.factor(C,hypers(2));  % factor analysis
        extras.loading_matrix = L;
        extras.indep_vars = psi(:);
        C = L*L' + diag(psi);
        
    case 'glasso'
        assert(length(hypers)==1)
        cove.set('max_latent',0)   % prevent latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),10,cove.set);
        extras.S = extras.S/scale;
        C = inv(extras.S);
        
    case 'L1'        
        assert(length(hypers)==1)
        C = L1precisionBCD(C/mean(diag(C)),hypers(1));
        
    case 'lv-glasso'
        assert(length(hypers)==2)
        cove.set('max_latent',inf)   % allow any number of latent variables
        scale = mean(diag(C));
        extras = cove.lvglasso(C/scale,hypers(1),hypers(2),cove.set);
        extras.S = extras.S/scale;  % scale back
        extras.L = extras.L/scale;
        [H,D] = svds(double(extras.L),sum(~~extras.eigL));
        extras.H = H*sqrt(D);
        C = inv(extras.S - extras.H*extras.H');        

    otherwise
        error 'unknown covariance estimator'
end

C = (C+C')/2;
