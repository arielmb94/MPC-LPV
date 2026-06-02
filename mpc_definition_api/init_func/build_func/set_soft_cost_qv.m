function mpc = set_soft_cost_qv(mpc)

if mpc.has_s_cnstr
    mpc.s_cnstr = init_gradSlack(mpc,mpc.s_cnstr);
end

if mpc.has_y_cnstr
    mpc.y_cnstr = init_gradSlack(mpc,mpc.y_cnstr);
end

if mpc.has_h_cnstr
    mpc.h_cnstr = init_gradSlack(mpc,mpc.h_cnstr);
end

end


function cnstr = init_gradSlack(mpc,cnstr)

if cnstr.min_limit
    if ~any(cnstr.qv_min)
        if ~isempty(mpc.Qe)
            cnstr.qv_min(:) = max(eig(mpc.Qe))*10;
        else
            cnstr.qv_min(:) = mpc.qv;
        end
    end
end

if cnstr.max_limit
    if ~any(cnstr.qv_max)
        if ~isempty(mpc.Qe)
            cnstr.qv_max(:) = max(eig(mpc.Qe))*10;
        else
            cnstr.qv_max(:) = mpc.qv;
        end
    end
end
end