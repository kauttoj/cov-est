function [mat,ind] = concat_cell(cell_array,upper)
%CONCAT_CELL Summary of this function goes here
%   Detailed explanation goes here

N = size(cell_array{1},1);
ind = nan;

if upper
    ind = find(triu(ones(N,N),1));
    mat = nan(length(ind),length(cell_array));
    for i=1:length(cell_array)
        mat(:,i)=cell_array{i}(ind);
    end    
else
    mat = cell_array{1};
    for i=2:length(cell_array)
        mat = cat(3,mat,cell_array{i});
    end
end

end

