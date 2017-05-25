function [P_u] = pp_step2(p)
    A = p.A; B = p.B; Q = p.Q;
    R = p.R; n = p.n;
    T = p.T; U = p.U;

    H = [A, (-B*(R\B') + T); (-Q - U), -A'];
    % [V, J] = jordan(H);
    % it is far more efficient to make use
    % of the generalized eigenvectors insted of
    % finding the jordan form because of the
    % slow implementation of MATLAB

    [V_, D] = eig(H);
    % need to swap the columns of the
    % eigenvectors to get a diagonal
    % matrix with the negative eigenvalues
    % followed by the positive eigenvalues
    V = p.sort_array(V_, D, n);

    W_0 = V(1 : n, n+1: end);
    Z_0 = V(n + 1 : end, n + 1 : end);
    P_u = Z_0/W_0;
end
