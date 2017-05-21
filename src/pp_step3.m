function [P] = pp_step3(p)
    B = p.B; R = p.R;
    n = p.n; T = p.T;
    A_b = p.A + p.a*eye(n);
    Q_b = p.Q_0 + p.U + p.DQ;

    setlmis([]);
    M_b = lmivar(1, [n 1]);
    W_b = lmivar(1, [n 1]);

    %  1st lmi, -1 means LMI < 0
    %  (1,1) block, 's' means -A_b*W_b - W_b*A_b^T
    lmiterm([-1 1 1 W_b], -A_b, 1, 's');
    %  (1,1) block, B*R^-1*B^T
    lmiterm([-1 1 1 0], B*(R\B'));
    %  (1,1) block, -T
    lmiterm([-1 1 1 0], -T);
    %  (2,2) block, Q_b^-1
    lmiterm([-1 2 2 0], Q_b);
    % (2,1) block, W_b
    lmiterm([-1 2 1 W_b], 1, 1);

    % 2nd lmi, -2 means LMI < 0
    % (1,1) block, M_b
    lmiterm([-2 1 1 M_b], 1, 1);
    % (2,1) block, I_n
    lmiterm([-2 2 1 0], eye(n));
    % (2,2) block, W_b
    lmiterm([-2 2 2 W_b], 1, 1);

    % complete the LMI framework setting
    lmis = getlmis;
    % minimize the trace of M
    c = mat2dec(lmis, eye(n), zeros(n));
    % set relative accuracy 1e-5
    options = [1e-5,0,0,0,0];
    % solve min problem
    [copt xopt] = mincx(lmis, c, options);
    % solution matrix M_b
    M_b = dec2mat(lmis, xopt, M_b);
    % solution matrix W_b
    W_b = dec2mat(lmis, xopt, W_b);
    P = inv(W_b);
end
