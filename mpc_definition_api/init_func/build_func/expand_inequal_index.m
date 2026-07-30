function [S_gi,S_vi] = expand_inequal_index(cnstr,S_gi,S_vi,has_vi,k)


if cnstr.min_limit
    S_gi = [S_gi;cnstr.g_min_index_k(:,k)];
    if has_vi
        S_vi = [S_vi;cnstr.v_min_index_k(:,k)];
    else
        zero_vec = zeros(length(cnstr.g_min_index_k(:,k)),1);
        S_vi = [S_vi;zero_vec];
    end
end

if cnstr.max_limit
    S_gi = [S_gi;cnstr.g_max_index_k(:,k)];
    if has_vi
        S_vi = [S_vi;cnstr.v_max_index_k(:,k)];
    else
        zero_vec = zeros(length(cnstr.g_max_index_k(:,k)),1);
        S_vi = [S_vi;zero_vec];
    end
end

end