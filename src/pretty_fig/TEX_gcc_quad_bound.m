function TEX_gcc_quad_bound()
    addpath('../');
    addpath('~/downloads/matlab2tikz-master/src/');
    [A, B, R, R_t_inv, P, Pnot, Gc] = gcc_quad_bound();
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
