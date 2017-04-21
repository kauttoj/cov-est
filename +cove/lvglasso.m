function out = lvglasso(SigmaO,alpha,beta,opts)
%{ 
Adapation of the PGADM algorithm as implemented in ADMM_R.m. 
Additional features:
1. Upon convergence, non-zero parameters are refitted unconstrainted
2. Small parameters are retested for significance using the AICc. 

"Alternating Direction Methods for Latent Variable Gaussian Graphical Model Selection", 
appearing in Neural Computation, 2013, by Ma, Xue and Zou, for solving
Latent Variable Gaussian Graphical Model Selection
min <R,SigmaO> - logdet(R) + alpha ||S||_1 + beta Tr(L)
s.t. R = S - L,  R positive definte,  L positive semidefinite

Adaptation: Dimitri Yatsenko,  June 2013

Original Authors: Shiqian Ma, Lingzhou Xue and Hui Zou
Copyright (C) 2013 Shiqian Ma, The Chinese University of Hong Kong
http://www1.se.cuhk.edu.hk/~sqma/ADMM-LVGLasso
%} 

% default options 
def.continuation = true;
def.verbose = 0;
def.num_continuation = 6;
def.eta = 0.25;
def.mu = size(SigmaO,1);  % rule of thumb from original paper
def.muf = 1e-2;
def.maxiter = 1e4;
def.stoptol = 1e-5;
def.tau = 0.6;
def.refit = false;
def.max_latent = inf;

% parameter opts overrides defaults
if nargin<4
    opts = struct;
end
for f = setdiff(fieldnames(def),fieldnames(opts))'
    opts.(f{1}) = def.(f{1});
end

% initialization
n = size(SigmaO,1);
s.S = eye(n); s.L = zeros(n,n); s.Lambda = zeros(n);
mu = opts.mu; eta = opts.eta; tau = opts.tau;
~opts.verbose || fprintf('ADMM: alpha=%1.3f  beta=%1.3f\n',alpha,beta); %#ok<VUNUS>
max_latent = min(n-1,opts.max_latent);

% Phase 1: L1-regularized sparsification (equavelent to original algorithm).
for iter = 1:opts.maxiter
    s = update(s,SigmaO,alpha,beta,mu,tau,max_latent);
    resid = norm(s.R-s.S+s.L,'fro')/max([1,norm(s.R,'fro'),norm(s.S,'fro'),norm(s.L,'fro')]);
    if resid < opts.stoptol || iter >= opts.maxiter, break, end
    if opts.continuation && mod(iter,opts.num_continuation)==0; mu = max(mu*eta,opts.muf); end;
end

if opts.refit
    % Phase 2: Release L1 constraint on non-zero parameters
    alpha = alpha.*ones(n);
    s.Lambda = zeros(n);
    for iter = iter:opts.maxiter
        alpha(~~abs(s.S)) = 0;  % release constraint on nonzero parameters
        s = update(s,SigmaO,alpha,beta,mu,tau,max_latent);
        resid = norm(s.R-s.S+s.L,'fro')/max([1,norm(s.R,'fro'),norm(s.S,'fro'),norm(s.L,'fro')]);
        if resid < opts.stoptol || iter >= opts.maxiter, break, end
        if opts.continuation && mod(iter,opts.num_continuation)==0; mu = max(mu*eta,opts.muf); end;
    end
end

% return results
resid = norm(s.R-s.S+s.L,'fro')/max([1,norm(s.R,'fro'),norm(s.S,'fro'),norm(s.L,'fro')]);
obj = objective(s,SigmaO,alpha,beta);  
out = s; out.obj = obj; out.resid = resid; out.iter = iter;



function s = update(s,SigmaO,alpha,beta,mu,tau,max_latent)
% the core of the AMDD algorithm

% update s.R
B = mu*SigmaO - mu*s.Lambda - s.S + s.L;
try
[U,D] = eig(B); 
catch
    % this is a hack that fixed a few cases of "eig did not not converge"
    [U,D] = eigs(B,size(B,1)-1);
end
d = diag(D);
s.eigR = (-d + sqrt(d.^2+4*mu))/2;
s.R = sym(U*diag(s.eigR)*U');

% update s.S and s.L
Gradpartial = s.S - s.L - s.R + mu*s.Lambda;
G = s.S - tau * Gradpartial; H = s.L + tau*Gradpartial;
s.S = sym(soft_shrink(G, tau*mu*alpha));
[U,D] = eig_(H,max_latent); s.eigL = max(0,diag(D)-tau*mu*beta);
s.L = sym(U*diag(s.eigL)*U');

% update s.Lambda
s.Lambda = sym(s.Lambda - (s.R-s.S+s.L)/mu);


function A = sym(A)
A = (A+A')/2;

function x = soft_shrink(z,tau)
x = sign(z).*max(abs(z)-tau,0);


function obj = objective(s,SigmaO,alpha,beta)
obj = sum(sum((s.S-s.L).*SigmaO)) - logDet(s.S-s.L) + sum(sum(alpha.*abs(s.S))) + beta*sum(s.eigL);


function [U,d] = eig_(H,max_latent)
% this is faster than eigs
try
    [U,d] = eig(H);
catch
    [U,d] = eigs(H,max_latent);
end
d = diag(d);
[~,ix] = sort(d,'descend');
ix = ix(1:max_latent);
U = U(:,ix);
d = diag(d(ix));

function [ld,A] = logDet(A)
d = eig(A);
if ~all(d>=0)
    ld = nan;
else
    ld = sum(log(d));
end











% 
% 
% function out = ADMM_R(SigmaO,alpha,beta,opts) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This file implements the PGADM algorithm described in 
% % "Alternating Direction Methods for Latent Variable Gaussian Graphical 
% % Model Selection", appeared in Neural Computation, 2013,
% % by Ma, Xue and Zou, for solving 
% % Latent Variable Gaussian Graphical Model Selection  
% % min <R,SigmaO> - logdet(R) + alpha ||S||_1 + beta Tr(L) 
% % s.t. R = S - L,  R positive definte,  L positive semidefinite 
% % 
% % Authors: Shiqian Ma, Lingzhou Xue and Hui Zou 
% % Copyright (C) 2013 Shiqian Ma, The Chinese University of Hong Kong  
% % Date: Jan 25, 2013
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % initialization    alpha = 0.05; beta = 0.25
% QUIET    = 1;
% ABSTOL   = 1e-5;
% RELTOL   = 1e-5;
% n = size(SigmaO,1); 
% R = eye(n,n); S = R; L = zeros(n,n); Lambda = zeros(n,n); 
% mu = opts.mu; eta = opts.eta; tau = opts.tau; 
% 
% for iter = 1: opts.maxiter
%     % update R
%     B = mu*SigmaO - mu*Lambda - S + L; 
%     [U,D] = mexeig(B); d = diag(D); 
%     eigR = (-d + sqrt(d.^2+4*mu))/2; 
%     R = U*diag(eigR)*U'; R = (R+R')/2; 
%     
%     S_old = S; 
%     L_old = L;
%     
%     % update S and L 
%     Gradpartial = S - L - R + mu*Lambda; 
%     G = S - tau * Gradpartial; H = L + tau * Gradpartial; 
%     S = soft_shrink(G, tau*mu*alpha); S = (S+S')/2; 
%     [U,D] = mexeig(H); d = diag(D)-tau*mu*beta; eigL = max(d,0); 
%     L = U*diag(eigL)*U'; L = (L+L')/2; 
%     
%     % update Lambda 
%     Lambda = Lambda - (R-S+L)/mu; Lambda = (Lambda + Lambda')/2; 
%     % diagnostics, reporting, termination checks
%     k = iter; 
%     history.objval(k)  = objective(R,SigmaO,eigR,S,eigL,alpha,beta); 
% 
%     history.r_norm(k)  = norm(R-S+L,'fro');
%     history.s_norm(k)  = norm(-([S;L] - [S_old;L_old])/mu);
% 
%     history.eps_pri(k) = sqrt(3*n*n)*ABSTOL + RELTOL*max(norm(R,'fro'), norm([S;-L],'fro'));
%     history.eps_dual(k)= sqrt(3*n*n)*ABSTOL + RELTOL*norm(Lambda,'fro');
% 
%     if ~QUIET
%         fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
%             history.r_norm(k), history.eps_pri(k), ...
%             history.s_norm(k), history.eps_dual(k), history.objval(k));
%     end
%  
%     % print stats
%     resid = norm(R-S+L,'fro')/max([1,norm(R,'fro'),norm(S,'fro'),norm(L,'fro')]); 
%     obj = history.objval(k);  
%     % check stop 
%     if resid < opts.stoptol 
%         out.R = R; out.S = S; out.L = L; out.obj = obj; out.eigR = eigR; out.eigL = eigL; out.resid = resid; out.iter = iter; return; 
%     end
%     if opts.continuation && mod(iter,opts.num_continuation)==0; mu = max(mu*eta,opts.muf); end;  
% end
% out.R = R; out.S = S; out.L = L; out.obj = obj; out.eigR = eigR; out.eigL = eigL; out.resid = resid; out.iter = iter;
% 
% 
% function x = soft_shrink(z,tau)
% x = sign(z).*max(abs(z)-tau,0); 
% 
% function obj = objective(R,SigmaO,eigR,S,eigL,alpha,beta)
% obj = sum(sum(R.*SigmaO)) - sum(log(eigR)) + alpha*sum(abs(S(:))) + beta*sum(eigL);
