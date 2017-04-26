function [ output_args ] = scikit_learn_lasso(X)
%SCIKIT_LEARN_LASSO Summary of this function goes here
%   Detailed explanation goes here


end


function cost = scikit_learn_cost(mle,precision_,alpha)
    p = size(precision_,1);    % neurons
    cost = - 2. * log_likelihood(mle, precision_) + p*log(2*pi);
    cost = cost + alpha*(sum(vectorize(abs(precision_))) - abs(trace(precision_)));
end

function log_likelihood_ = log_likelihood(emp_cov, precision)
    p = size(precision,1);
    log_likelihood_ = - sum(vectorize(emp_cov*precision)) + cove.logDet(precision); % same as above
    log_likelihood_ = log_likelihood_ - p*log(2*pi);
    log_likelihood_ = log_likelihood_/2;
end