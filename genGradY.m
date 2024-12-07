function [gradYmin,gradYmax] = genGradY(C,D,N,Nx,Nu,Ny,nx,nu,ny)
%optimization variables are stored as [u0 x1 u1 x2 ... xn-1 un-1 xn]^T

gradYmax = zeros(Nx+Nu,Ny);
for k = 0:N-1

    switch k

        case 0
            if D == 0
                % Do Nothing
            else
                % gradYmax = D'
                gradYmax(1 : nu,1:ny) = D';
            end

        otherwise

            if D == 0 % -> Ny is (N-1)*ny
                % gradYmax = [C';D']
                gradYmax(nu + 1 + nx*(k-1) + nu*(k-1): nu + nx*(k)+ nu*(k),...
                    (ny*(k)+1:ny*(k+1))-ny) = [ C' ; D'];
            else      % -> Ny is N*ny
                % gradYmax = [C';D']
                gradYmax(nu + 1 + nx*(k-1) + nu*(k-1): nu + nx*(k)+ nu*(k),...
                    ny*(k)+1:ny*(k+1)) = [ C' ; D'];
            end
    end

end

gradYmin = -gradYmax;
end