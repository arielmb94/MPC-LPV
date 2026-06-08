function v_full = fill_vec(v_full, v, fill_style)

[~, step_v] = size(v);
[~, step_full] = size(v_full);

% Copy available stages
v_full(:,1:step_v) = v;

% Fill remaining stages
if fill_style == 0
    % Zero fill
    v_full(:,step_v+1:step_full) = 0;
else
    % Repeat last stage
    for i = step_v+1:step_full
        v_full(:,i) = v(:,step_v);
    end
end

end