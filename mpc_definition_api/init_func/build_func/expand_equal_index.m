function [start_index,index_k_mat] = expand_equal_index(dim,start_index,index_k_mat)

index_k_mat = [index_k_mat;start_index:start_index+dim-1];
start_index = start_index+dim;
end