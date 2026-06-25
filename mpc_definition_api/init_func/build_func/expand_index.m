function [start_index,index,index_k] = expand_index(dim,start_index,i,index)

index_k = [start_index:start_index+dim(i)-1]';
index = [index; index_k];

start_index = start_index+dim(i);
end