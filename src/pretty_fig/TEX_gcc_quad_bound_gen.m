function TEX_gcc_quad_bound_gen()
    addpath('../');
    addpath('~/downloads/matlab2tikz-master/src/');
    [A, B, R, P, Pnot, Gc, e] = gcc_quad_bound_gen();
    makeplot(A, B, R, P, Pnot, Gc, e);
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
