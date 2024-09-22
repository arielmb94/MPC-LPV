function hess = hess_linear_Ind(fi,grad_fi)

    % Hess(Phi) = sum((grad(fi)*grad(fi)^T)/fi^2) 
  
    [n,m] = size(grad_fi);

    hess = zeros(n);
    for i = 1:m

        hess = hess + (grad_fi(:,i)*grad_fi(:,i)')/(fi(i)^2);

    end

end