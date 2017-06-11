function [S] = pp_step1(p)
    A = p.A; B = p.B;
    Q_0 = p.Q_0; R = p.R; n = p.n;
    k = p.k; D = p.D; E = p.E;

    setlmis([]);
    M = lmivar(1, [n 1]);
    W = lmivar(1, [n 1]);
    %  k = 2
    S = lmivar(1, repmat([1 0], k, 1));

    %  1st lmi, -1 means LMI < 0
    %  (1,1) block, 's' means -A*W - W*A^T
    lmiterm([-1 1 1 W], -A, 1, 's');
    %  (1,1) block, B*R^-1*B^T
    lmiterm([-1 1 1 0], B*(R\B'));
    %  (1,1) block, D*S*D^T
    lmiterm([-1 1 1 S], -D, D');
    %  (2,2) block, S
    lmiterm([-1 2 2 S], 1, 1);
    %  (3,3) block, Q_0^-1
    lmiterm([-1 3 3 0], Q_0);
    % (2,1) block, E*W
    lmiterm([-1 2 1 W], E, 1);
    % (3,1) block, W
    lmiterm([-1 3 1 W], 1, 1);

    % 2nd lmi, -2 means LMI < 0
    % (1,1) block, M
    lmiterm([-2 1 1 M], 1, 1);
    % (2,1) block, I_n
    lmiterm([-2 2 1 0], eye(n));
    % (2,2) block, W
    lmiterm([-2 2 2 W], 1, 1);

    % complete the LMI framework setting
    lmis = getlmis;
    % minimize the trace of M
    c = mat2dec(lmis, eye(n), zeros(n), zeros(k));
    % set relative accuracy 1e-5
    options = [1e-5,0,0,0,0];
    % solve min problem
    [copt, xopt] = mincx(lmis, c, options);
    % solution matrix M
    M = dec2mat(lmis, xopt, M);
    % solution matrix W
    W = dec2mat(lmis, xopt, W);
    % solution matrix S
    S = dec2mat(lmis, xopt, S);
    P = inv(W);
end
