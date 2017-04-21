function res = edgeoverlap(mat1,mat2,percentage)
%EDGEOVERLAP Summary of this function goes here
%   Detailed explanation goes here

N_cells = size(mat1,1);

IND = find(triu(ones(N_cells,N_cells),1));

mat1 = nodiag(mat1);
mat2 = nodiag(mat2);

th1 = prctile(abs(mat1(IND)),100-percentage);
th2 = prctile(abs(mat2(IND)),100-percentage);

net1 = double(abs(mat1)>th1);
net2 = double(abs(mat2)>th2);

res = sum(sum(net1.*net2))/(0.5*sum(sum(net1))+0.5*sum(sum(net2)));

end

