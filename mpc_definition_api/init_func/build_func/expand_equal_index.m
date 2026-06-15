function [start_index,index_k_mat,index_k_type] = expand_equal_index(dim,...
                                                    start_index,index_k_mat,...
                                                    index_k_type)
index_k = (start_index:start_index+dim-1)';
index_k_mat = [index_k_mat index_k];
index_k_type = [index_k_type;index_k];
start_index = start_index+dim;
end