function M_full = fill_mat(M_full, M, fill_style)

step_m = size(M, 3);
step_full = size(M_full, 3);

% Copy available stages 
M_full(:,:, 1:step_m) = M;

% Fill remaining stages
if step_m < step_full
    if fill_style == 0
        % Zero fill remaining stages
        M_full(:,:, step_m+1:step_full) = 0;
    else
        % Repeat last available stage
        for k = step_m+1:step_full
            M_full(:,:, k) = M(:,:, step_m);
        end
    end
end

end