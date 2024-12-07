% Computes gradient and hessian for cost terms of the form:
% 0.5*y'Qy
% where y is a linear output equation:
% y = Cs + Du + Dd
%
% with:
%
% grad = Grad(y)*Q*y
% hess = Grad(y)*Q*Delta_y'
%
% For error cost:
% e = r-y = r-(Cs + Du + Dd)
% grad = -[C'; D']*Q*e
% hess = [C' ; D']*Q*[C D]
%
% For general linear perf. cost:
% z = Cs + Du + Dd
% grad = [C'; D']*Q*z
% hess = [C' ; D']*Q*[C D]

function [gradQ,hess] = genLinOutGradHess(Q,C,D,N,Nx,Nu,Ny,nx,nu,ny)

% grad = Grad(y)*Q*y
% gradQ = Grad(y)*Q
gradQ = zeros(Nx+Nu,Ny);
hess = zeros(Nx+Nu,Nx+Nu);

for k = 0:N-1

    switch k

        case 0
            
            if D==0
                % Do nothing
            else
                % gradQ = D'*Q
                gradQ(1:nu, 1:ny) = D'*Q;           
                % hess = gradQ*D
                hess(1:nu, 1:nu) = gradQ(1:nu, 1:ny)*D;
            end

        otherwise

            if D==0 % -> Ny is (N-1)*ny
                
                % gradQ = [C';D']*Q
                gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    (ny*(k)+1:ny*(k+1))-ny) = ...
                    [ C' ; D']*Q;
                
                % hess = gradQ*[C D]
                hess(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k))...
                    = gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    (ny*(k)+1:ny*(k+1))-ny)*[C D];
            
            else    % -> Ny is N*ny

                % gradQ = [C';D']*Q
                gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    ny*(k)+1:ny*(k+1)) = ...
                    [ C' ; D']*Q;

                % hess = gradQ*[C D]
                hess(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k))...
                    = gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                    ny*(k)+1:ny*(k+1))*[C D];
            end

    end
    
end

end