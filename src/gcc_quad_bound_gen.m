function [A, B, R, P, Pnot, Gc, e] = gcc_quad_bound_gen();
    A = [-1, 2; -1, -1];
    B = [0; 1];
    Q = eye(length(A));
    R = 10;

    d_t1 = [0.9; 0.3];
    e_t1 = [0.6; 0.2];
    d_t2 = [0.1; 0.3];
    e_t2 = [0.4; 0.9];
    f_t1 = [0.3; 0.2];
    g_t1 = 0.6;

    DA = @(r1, r2) r1*d_t1*e_t1' + r2*d_t2*e_t2';
    DB = @(p1) p1*f_t1*g_t1';

    T = d_t1*d_t1' + d_t2*d_t2';
    U = e_t1*e_t1' + e_t2*e_t2';
    V = g_t1*g_t1';
    W = f_t1*f_t1';

    e = 1;
    M = @(e) (1/e)*B*(R\B') ...
	- (1/e^2)*B*(R\(V*(R\B'))) - W - T;
    H = [A, -M(e); -Q - U, -A'];
    lambdaH = eig(H);

    Pnot = care(A, B, Q, R);
    [P, fval] = fsolve(@(p) gare(p, A, M(e), Q, U), Pnot);
    if all(eig(P) > eps)
	Gc = @(r1, r2, p1) A + DA(r1, r2) ...
	    - (1/e)*(B + DB(p1))*(R\(B'*P));
	makeplot(A, B, R, P, Pnot, Gc, e);
    end
end

function f = gare(p, A, M, Q, U)
    P = reshape(p(1:length(A)^2), size(A));
    F = P*A + A'*P - P*M*P + Q + U;
    f = F(:);
end

function dxdt = state(~, x, Gc)
    dxdt = Gc*x;
end

function makeplot(A, B, R, P, Pnot, Gc, e)
    u = @(P) -R\(B'*P);
    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
	state(t, x, A - B*(R\(B'*Pnot))), [0 7], [3; -2]);
    plot(t, x, t, u(Pnot)*x','LineWidth',1.0);
    title('Initial System')
    legend('x_1(t)', 'x_2(t)', 'u(t)');
    grid on;

    figure(2)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(1, 1, 1)), [0 10], [3; -2]);
    plot(t, x, t, (1/e)*u(P)*x','LineWidth',1.0);
    title('Uncertain System r1 = 1, r2 = 1, p1 = 1')
    legend('x_1(t)', 'x_2(t)', 'u(t)');
    grid on;

    figure(3)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(-1, -1, -1)), [0 5], [3; -2]);
    plot(t, x, t, (1/e)*u(P)*x','LineWidth',1.0);
    title('Uncertain System r1 = -1, r2 = -1, p1 = -1')
    legend('x_1(t)', 'x_2(t)', 'u(t)');
    grid on;
end
