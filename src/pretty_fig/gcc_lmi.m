function [P] = gcc_lmi()
    A = [-1, 2; -1, -1];
    B = [0; 1];
    Q = eye(length(A));
    Q_0 = inv(Q);
    R = 10;

    d_t1 = [0.9; 0.3];
    e_t1 = [0.6; 0.2];
    d_t2 = [0.1; 0.3];
    e_t2 = [0.4; 0.9];
    f_t1 = [0.3; 0.2];
    g_t1 = 0.6;

    n = length(A);
    k = 2;
    l = 1;

    D = [d_t1, d_t2];
    E = [e_t1, e_t2];
    F = f_t1;
    G = g_t1;

    setlmis([]);
    M_t = lmivar(1, [n 1]);
    P_t = lmivar(1, [n 1]);
    %  k = 2
    S_t = lmivar(1, [1 0; 1 0]);
    %  l = 1
    T_t = lmivar(1, [1 0]);
    d = lmivar(1, [1 0]);

    %  1st lmi, -1 means LMI < 0
    %  (1,1) block, 's' means -A*P_t - P_t*A^T
    lmiterm([-1 1 1 P_t], -A, 1, 's');
    %  (1,1) block, d*B*R^-1*B^T
    lmiterm([-1 1 1 d], 1, B*(R\B'));
    %  (1,1) block, D*S_t*D^T
    lmiterm([-1 1 1 S_t], -D, D');
    %  (1,1) block, F*T_t*F^T
    lmiterm([-1 1 1 T_t], -F, F');
    %  (2,2) block, S_t
    lmiterm([-1 2 2 S_t], 1, 1);
    %  (3,3) block, T_t
    lmiterm([-1 3 3 T_t], 1, 1);
    %  (4,4) block, Q_0^-1
    lmiterm([-1 4 4 0], Q_0);
    % (2,1) block, E*P_t
    lmiterm([-1 2 1 P_t], E, 1);
    % (3,1) block, d*G*R^-1*B^T
    lmiterm([-1 3 1 d], 1, G*(R*B'));
    % (4,1) block, P_t
    lmiterm([-1 4 1 P_t], 1, 1);

    % 2nd lmi, -2 means LMI < 0
    % (1,1) block, M_t
    lmiterm([-2 1 1 M_t], 1, 1);
    % (2,1) block, I_n
    lmiterm([-2 2 1 0], eye(n));
    % (2,2) block, P_t
    lmiterm([-2 2 2 P_t], 1, 1);

    % complete the LMI framework setting
    lmis = getlmis;
    % minimize the trace of M_t
    c = mat2dec(lmis, eye(n), zeros(n), zeros(k), zeros(l), 0);
    % set relative accuracy 1e-5
    options = [1e-5,0,0,0,0];
    % solve min problem
    [copt xopt] = mincx(lmis, c, options);
    % solution matrix M_t
    M_t = dec2mat(lmis, xopt, M_t);
    % solution matrix P_t
    P_t = dec2mat(lmis, xopt, P_t);
    % solution matrix S_t
    S_t = dec2mat(lmis, xopt, S_t);
    % solution matrix T_t
    T_t = dec2mat(lmis, xopt, T_t);
    % solution matrix d
    d = dec2mat(lmis, xopt, d);
    P = inv(P_t);

    DA = @(r1, r2) ...
	r1*S_t(1,1)^(0.5)*d_t1*S_t(1,1)^(-0.5)*e_t1' ...
	+ r2*S_t(2,2)^(0.5)*d_t2*S_t(2,2)^(-0.5)*e_t2';
    DB = @(p1) ...
	p1*T_t(1,1)^(0.5)*f_t1*T_t(1,1)^(-0.5)*g_t1';

    Gc = @(r1, r2, p1) ...
	A + DA(r1, r2) - (B + DB(p1))*(R\(B'*P));
    Pnot = care(A, B, Q, R);
    makeplot(A, B, R, P, Pnot, Gc, d);
end

function dxdt = state(~, x, Gc)
    dxdt = Gc*x;
end

function makeplot(A, B, R, P, Pnot, Gc, d)
    u = @(P) -R\(B'*P);
    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
	state(t, x, A - B*(R\(B'*Pnot))), [0 7], [3; -2]);
    plot(t, x, t, u(Pnot)*x','LineWidth',1.0);
    makepretty('gcc_lmi1')

    figure(2)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(1, 1, 1)), [0 10], [3; -2]);
    plot(t, x, t, d*u(P)*x','LineWidth',1.0);
    makepretty('gcc_lmi2')

    figure(3)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(-1, -1, -1)), [0 5], [3; -2]);
    plot(t, x, t, d*u(P)*x','LineWidth',1.0);
    makepretty('gcc_lmi3')
end

function makepretty(string)
    xlabel('Time (seconds)','interpreter','latex');
    ylabel('Amplitude','interpreter','latex');
    grid on;
    l = legend(['$x_1(t)$'], ['$x_2(t)$'], ['$u(t)$']);
    set(l,'Interpreter','Latex');
    cleanfigure;
    matlab2tikz(strcat(string, '.tex'));
end
