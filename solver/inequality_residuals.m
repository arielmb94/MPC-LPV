function [ri_0,ri_k,ri_ter] = inequality_residuals(mpc,N,ri_0,ri_k,ri_ter,...
                    g_0,g,g_ter,v_0,v,v_ter,...
                    ng_k,nv_k,v_rows_0,v_rows_k,...
                    s_cnstr,u_cnstr,du_cnstr,y_cnstr,h_cnstr)

%%
% Inequality residuals are assembled directly from materialized signals.
% A lower-bound row is q_min - q + g - E*v and an upper-bound row is
% q - q_max + g - E*v.

% k = 0
if ng_k(1)
    ri_0(:) = g_0;
    if nv_k(1)
        ri_0(v_rows_0) = ri_0(v_rows_0) - v_0;
    end
end

if ng_k(2)
    ri_k(:,:) = g;
    if nv_k(2)
        ri_k(v_rows_k,:) = ri_k(v_rows_k,:) - v;
    end
end

if ng_k(3)
    ri_ter(:) = g_ter- v_ter;
end

if mpc.has_s_cnstr
    [ri_k,ri_ter] = s_residual(ri_k,ri_ter,s_cnstr.min_limit,s_cnstr.max_limit,...
                            s_cnstr.min,s_cnstr.max,...
                            s_cnstr.min_ineqRow_k,s_cnstr.min_ineqRow_ter,...
                            s_cnstr.max_ineqRow_k,s_cnstr.max_ineqRow_ter,...
                            mpc.s,N);
end

if mpc.has_u_cnstr
    [ri_0,ri_k] = u_residual(ri_0,ri_k,u_cnstr.min_limit,u_cnstr.max_limit,...
                            u_cnstr.min,u_cnstr.max,...
                            u_cnstr.min_ineqRow_0,u_cnstr.min_ineqRow_k,...
                            u_cnstr.max_ineqRow_0,u_cnstr.max_ineqRow_k,...
                            mpc.u,N);
end

if mpc.has_du_cnstr
    [ri_0,ri_k] = u_residual(ri_0,ri_k,du_cnstr.min_limit,du_cnstr.max_limit,...
                            du_cnstr.min,du_cnstr.max,...
                            du_cnstr.min_ineqRow_0,du_cnstr.min_ineqRow_k,...
                            du_cnstr.max_ineqRow_0,du_cnstr.max_ineqRow_k,...
                            mpc.du,N);
end

if mpc.has_y_cnstr
    [ri_0,ri_k,ri_ter] = y_residual(ri_0,ri_k,ri_ter,...
                            y_cnstr.min_limit,y_cnstr.max_limit,...
                            y_cnstr.min_0,y_cnstr.min_k,y_cnstr.min_ter,...
                            y_cnstr.max_0,y_cnstr.max_k,y_cnstr.max_ter,...
                            y_cnstr.min_ineqRow_0,y_cnstr.min_ineqRow_k,y_cnstr.min_ineqRow_ter,...
                            y_cnstr.max_ineqRow_0,y_cnstr.max_ineqRow_k,y_cnstr.max_ineqRow_ter,...
                            y_cnstr.use_k0,y_cnstr.use_ter,...
                            mpc.y_0,mpc.y,mpc.y_ter);
end


if mpc.has_h_cnstr
    [ri_0,ri_k,ri_ter] = y_residual(ri_0,ri_k,ri_ter,...
                            h_cnstr.min_limit,h_cnstr.max_limit,...
                            h_cnstr.min_0,h_cnstr.min_k,h_cnstr.min_ter,...
                            h_cnstr.max_0,h_cnstr.max_k,h_cnstr.max_ter,...
                            h_cnstr.min_ineqRow_0,h_cnstr.min_ineqRow_k,h_cnstr.min_ineqRow_ter,...
                            h_cnstr.max_ineqRow_0,h_cnstr.max_ineqRow_k,h_cnstr.max_ineqRow_ter,...
                            h_cnstr.use_k0,h_cnstr.use_ter,...
                            mpc.h_0,mpc.h,mpc.h_ter);
end

end

function [ri_0,ri_k] = u_residual(ri_0,ri_k,min_limit,max_limit,min,max,...
                            min_ineqRow_0,min_ineqRow_k,...
                            max_ineqRow_0,max_ineqRow_k,u,N)

if min_limit
    row = min_ineqRow_0;
    
    ri_u = min - u;

    ri_0(row) = ri_0(row) + ri_u(:,1);

    row = min_ineqRow_k;
    ri_k(row,:) = ri_k(row,:) + ri_u(:,2:N);
end
if max_limit
    row = max_ineqRow_0;

    ri_u = u-max;

    ri_0(row) = ri_0(row) + ri_u(:,1);

    row = max_ineqRow_k;
    ri_k(row,:) = ri_k(row,:) + ri_u(:,2:N);
end

end

function [ri_k,ri_ter] = s_residual(ri_k,ri_ter,min_limit,max_limit,min,max,...
                            min_ineqRow_k,min_ineqRow_ter,...
                            max_ineqRow_k,max_ineqRow_ter,s,N)

if min_limit

    row = min_ineqRow_k;

    ri_s = min-s;

    ri_k(row,:) = ri_k(row,:) + ri_s(:,1:N-1);

    row = min_ineqRow_ter;
    ri_ter(row) = ri_ter(row) + ri_s(:,N);
end
if max_limit
    row = max_ineqRow_k;

    ri_s = s-max;

    ri_k(row,:) = ri_k(row,:) + ri_s(:,1:N-1);

    row = max_ineqRow_ter;
    ri_ter(row) = ri_ter(row) + ri_s(:,N);
end

end

function [ri_0,ri_k,ri_ter] = y_residual(ri_0,ri_k,ri_ter,...
                            min_limit,max_limit,...
                            min_0,min_k,min_ter,max_0,max_k,max_ter,...
                            min_ineqRow_0,min_ineqRow_k,min_ineqRow_ter,...
                            max_ineqRow_0,max_ineqRow_k,max_ineqRow_ter,...
                            use_k0,use_ter,y_0,y,y_ter)

if min_limit
    if use_k0
        row = min_ineqRow_0;
        ri_0(row) = ri_0(row) + min_0-y_0;
    end
    
    row = min_ineqRow_k;
    ri_k(row,:) = ri_k(row,:) + min_k-y;

    if use_ter
        row = min_ineqRow_ter;
        ri_ter(row) = ri_ter(row) + min_ter-y_ter;
    end
end
if max_limit
    if use_k0
        row = max_ineqRow_0;
        ri_0(row) = ri_0(row) + y_0-max_0;
    end
    
    row = max_ineqRow_k;
    ri_k(row,:) = ri_k(row,:) + y-max_k;

    if use_ter
        row = max_ineqRow_ter;
        ri_ter(row) = ri_ter(row) + y_ter-max_ter;
    end
end

end
