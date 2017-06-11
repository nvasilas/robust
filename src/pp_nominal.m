function [Gnom, Pnom, error_flag] = pp_nominal(p)
    A = p.A; B = p.B;
    Q = p.Q; R = p.R;
    n = p.n; a = p.a;

    H = [A, -B*(R\B'); -Q, -A'];
    % [V, J] = jordan(H);
    % it is far more efficient to make use
    % of the generalized eigenvectors insted of
    % finding the jordan form because of the
    % slow implementation of MATLAB

    [V_, D] = eig(H);
    if 2*n ~= length(unique(diag(D)))
        error_flag = true;
    else
        error_flag = false;
    end
    % need to swap the columns of the
    % eigenvectors to get a diagonal
    % matrix with the negative eigenvalues
    % followed by the positive eigenvalues
    V = p.sort_array(V_, D, n);

    W_0 = V(1 : n, n+1: end);
    Z_0 = V(n + 1 : end, n + 1 : end);
    P_u = Z_0/W_0;
    DQ = -2*a*P_u;

    Pnom = care(A + a*eye(n), B, Q + DQ, R);
    Gnom = A - B*(R\(B'*Pnom));
end
