function x = solveLinearSystemLDL(A, b)
    % "solveLinearSystemLDL" Solves Ax = b using LDL' factorization.
    %   A: Symmetric positive definite matrix
    %   b: Right-hand side vector
    %   x: Solution vector
    
    [L, D] = ldlFactorization(A);
    y = forwardSubstitution(L, b);
    z = diagonalSolve(D, y);
    x = backSubstitution(L', z);
end
function [L, D] = ldlFactorization(A)
    % "ldlFactorization" Performs LDL' decomposition.
    %   A: Symmetric positive definite matrix
    %   L: Lower triangular matrix with unit diagonal
    %   D: Diagonal matrix
    
    n = size(A, 1);
    L = eye(n);
    D = zeros(n);
    
    for i = 1:n
        D(i, i) = A(i, i) - L(i, 1:i-1) * (D(1:i-1, 1:i-1) * L(i, 1:i-1)');
        for j = i+1:n
            L(j, i) = (A(j, i) - L(j, 1:i-1) * (D(1:i-1, 1:i-1) * L(i, 1:i-1)')) / D(i, i);
        end
    end
end
function y = forwardSubstitution(L, b)
    % "forwardSubstitution" Solves Ly = b for y (L is lower triangular).
    n = length(b);
    y = zeros(n, 1);
    
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end
end
function z = diagonalSolve(D, y)
    % "diagonalSolve" Solves Dz = y for z (D is diagonal).
    z = y ./ diag(D);
end
function x = backSubstitution(U, y)
    % "backSubstitution" Solves Ux = y for x (U is upper triangular).
    n = length(y);
    x = zeros(n, 1);
    
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end