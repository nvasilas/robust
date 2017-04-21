function [P, lambdaH] = gcc_quad_bound_gen()
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
    makepretty('gcc_quad_bound_gen1')

    figure(2)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(1, 1, 1)), [0 10], [3; -2]);
    plot(t, x, t, (1/e)*u(P)*x','LineWidth',1.0);
    makepretty('gcc_quad_bound_gen2')

    figure(3)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(-1, -1, -1)), [0 5], [3; -2]);
    plot(t, x, t, (1/e)*u(P)*x','LineWidth',1.0);
    makepretty('gcc_quad_bound_gen3')
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
