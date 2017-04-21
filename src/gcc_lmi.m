function [P] = gcc_lmi()
    A = [-1, 2; -1, -1];
    B = [0; 1];
    Q = eye(length(A));
    R = 10;

    d_1 = [0.9; 0.3];
    e_1 = [0.6; 0.2];
    d_2 = [0.1; 0.3];
    e_2 = [0.4; 0.9];
    f_1 = [0.3; 0.2];
    g_1 = 0.6;

    n = length(A);
    k = 2;
    l = 1;

    D = [d_t1, d_t2];
    E = [e_t1, e_t2];
    F = f_t1;
    G = g_t1;

    setlmis([]);
    M_t = lmivar(1,[n 1]);
    P_t = lmivar(1,[n 1]);
    S_t = lmivar(1,[1 0; 1 0]); %k = 2
    T_t = lmivar(1,[1 0]); %l = 1
    d = lmivar(1,[1 0]);
end
