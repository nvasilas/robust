function [P] = gcc_quad_bound()
    A = [-1, 0; 1, -2];
    B = [0; 1];
    Q = eye(length(A));
    R = 10;

    A_t1 = [5.5, 0];
    A_t2 = [0, 2.0];
    B_t1 = 2.5;
    DA = @(r1, r2) B*(A_t1*r1 + A_t2*r2);
    DB = @(p1) B*B_t1*p1;
    Pnot = care(A, B, Q, R);

    a = 1;
    b = 1;

    theta_ = (1/a^2 + 1/a^2)*eye(size(B, 2));
    ksi_ = a^2*(A_t1')*A_t1 + a^2*(A_t2')*A_t2;
    phi_ = (1/b^2 + 1/b^2)*eye(size(B, 2));
    psi_ = b^2*R\(B_t1'*(B_t1/R));

    Q_t = Q + ksi_;
    R_t_inv = inv(R) - theta_ + phi_ + psi_;
    P = care(A, B, Q_t, inv(R_t_inv));
    Gc = @(r1, r2, p1) ...
	A + DA(r1, r2) - (B + DB(p1))*R_t_inv*B'*P;
    makeplot(A, B, R, R_t_inv, P, Pnot, Gc);
end

function dxdt = state(~, x, Gc)
    dxdt = Gc*x;
end

function makeplot(A, B, R, R_t_inv, P, Pnot, Gc)
    Knot = -R\(B'*Pnot);
    K = -R_t_inv*B'*P;
    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
	state(t, x, A - B*(R\(B'*Pnot))), [0 10], [-1; 2]);
    plot(t, x, t, Knot*x','LineWidth',1.0);
    makepretty('gcc_quad_bound1')

    figure(2)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(1, 1, 1)), [0 15], [-1; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    makepretty('gcc_quad_bound2')

    figure(3)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(-1, -1, -0.5)), [0 10], [-1; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    makepretty('gcc_quad_bound3')
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
