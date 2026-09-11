% INIT_DISCRETIZE_SYSTEM discretizes continuous-time state-space matrices.
    % 
    % Inputs:
    %   Ac, Bc, Bdc : Continuous-time system matrices
    %   Ts          : Sampling time
    %   method      : 'forward' (Forward Euler), 'backward' (Backward Euler), 
    %                 or 'tustin' (Bilinear transform). Default is 'forward'.
    % Outputs: 
    %   Ad, Bd, Bdd: Discrete-time system matrices
function [Ad, Bd, Bdd] = init_discretize_system(Ac, Bc, Bdc, Ts, method)
    
    if nargin < 5
        method = 'forward'; % Default method
    end
    
    if isempty(Bdc)
        Bdc = zeros(size(Ac,1), 1);
    end
    
    nx = size(Ac, 1);
    I = eye(nx);
    
    switch lower(method)
        case 'forward'
            % Forward Euler approximation
            Ad = I + Ts * Ac;
            Bd = Ts * Bc;
            Bdd = Ts * Bdc;
            
        case 'backward'
            % Backward Euler approximation
            inv_term = inv(I - Ts * Ac);
            Ad = inv_term;
            Bd = inv_term * (Ts * Bc);
            Bdd = inv_term * (Ts * Bdc);
            
        case 'tustin'
            % Tustin (Bilinear) approximation
            inv_term = inv(I - (Ts/2) * Ac);
            Ad = inv_term * (I + (Ts/2) * Ac);
            Bd = inv_term * (Ts * Bc);
            Bdd = inv_term * (Ts * Bdc);
            
        otherwise
            error('Unknown discretization method. Choose ''forward'', ''backward'', or ''tustin''.');
    end
end