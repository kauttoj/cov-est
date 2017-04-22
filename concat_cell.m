function mat = concat_cell(cell_array)
%CONCAT_CELL Summary of this function goes here
%   Detailed explanation goes here

mat = cell_array{1};
for i=2:length(cell_array)
    mat = cat(3,mat,cell_array{i});
end

end

