function [start_index,index,index_k_mat] = expand_index(dim,start_index,i,index,index_k_mat)
index_k = start_index:start_index+dim(i)-1;
index_k_mat = [index_k_mat;index_k];
index = [index index_k];
start_index = start_index+dim(i);
end