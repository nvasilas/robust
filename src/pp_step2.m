function [P_u] = pp_step2(p)
    A = p.A; B = p.B; Q = p.Q;
    R = p.R; n = p.n;
    T = p.T; U = p.U;

    H = [A, (-B*(R\B') + T); (-Q - U), -A'];
    [V, J] = jordan(H);
    W_0 = V(1 : n, n+1: end);
    Z_0 = V(n + 1 : end, n + 1 : end);
    P_u = Z_0/W_0;
end