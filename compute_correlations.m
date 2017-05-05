function [cormat,cormat_part] = compute_correlations(DATA)

% DATA = cell x file x stimulus x repeat x bin
N_cells = size(DATA,1);
N_file = size(DATA,2);
N_stim = size(DATA,3);
N_repeat = size(DATA,4);
N_bin = size(DATA,5);

DATA = permute(DATA,[1,3,5,2,4]);
DATA = reshape(DATA,N_cells,N_stim,N_bin,N_repeat*N_file);

M = mean(DATA,4);

DATA = bsxfun(@minus,DATA,M); % remove mean

DATA = permute(DATA,[1,2,4,3]);

cormat = 0;
cormat_part = 0;
for i = 1:N_stim
    dat = squeeze(DATA(:,i,:,:));
    dat = reshape(dat,N_cells,[]);
    
    c = corr(dat'); % correlations
    p = -corrcov(inv(cov(dat'))); % partial correlations
    
    assert(nnz(isnan(c))==0);    
    
    cormat = cormat + c/N_stim;
    cormat_part = cormat_part + p/N_stim;
end

end
