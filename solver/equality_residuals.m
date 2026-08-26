function [rp_0,rp_k] = equality_residuals(rp_0,rp_k,A,B,beq_0,beq_k,u,s,has_du,...
                                   su,s_idx,su_idx,N)


% k = 0
rp_0(s_idx) = B(:,:,1)*u(:,1)-s(:,1)-beq_0;

for k = 1:N-1
    % [A B -I][s u s+]'-beq
    rp_k(s_idx,k) = A(:,:,k+1)*s(:,k)+B(:,:,k+1)*u(:,k+1)-s(:,k+1)...
                            -beq_k(:,k);
end

if has_du
    % this is the control input delay equality condition: su+ = u_prev
    % [I -I][u su+]' = 0
    % k = 0
    rp_0(su_idx) = u(:,1)-su(:,1);
    % k = 1...N
    rp_k(su_idx,:) = u(:,2:N)-su(:,2:N);
end

end