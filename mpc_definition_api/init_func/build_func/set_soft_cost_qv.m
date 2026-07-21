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

if ~isempty(mpc.Qe)
    qv = max(eig(mpc.Qe))*10;
else
    qv = mpc.qv;
end

if cnstr.min_limit
    if cnstr.use_k0
        if ~any(cnstr.qv_min_0)
            cnstr.qv_min_0(:) = qv;
        end
    end
    if ~any(cnstr.qv_min)
        cnstr.qv_min(:) = qv;
    end
    if cnstr.use_ter
        if ~any(cnstr.qv_min_ter)
            cnstr.qv_min_ter(:) = qv;
        end
    end
end

if cnstr.max_limit
    if cnstr.use_k0
        if ~any(cnstr.qv_max_0)
            cnstr.qv_max_0(:) = qv;
        end
    end
    if ~any(cnstr.qv_max)
        cnstr.qv_max(:) = qv;
    end
    if cnstr.use_ter
        if ~any(cnstr.qv_max_ter)
            cnstr.qv_max_ter(:) = qv;
        end
    end
end

end