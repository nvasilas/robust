function TEX_gcc_linear_bound()
    addpath('../');
    addpath('~/downloads/matlab2tikz-master/src/');
    [A, B, R, P, Pnot, Gc] = gcc_linear_bound();
    makeplot(A, B, R, P, Pnot, Gc);
end

function dxdt = state(~, x, Gc)
    dxdt = Gc*x;
end

function makeplot(A, B, R, P, Pnot, Gc)
    u = @(P) -R\(B'*P);
    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
       	state(t, x, A - B*(R\(B'*Pnot))), [0 10], [1; 2]);
    plot(t, x, t, u(Pnot)*x','LineWidth',1.0);
    makepretty('gcc_linear_bound1')

    figure(2)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(1, 1)), [0 15], [1; 2]);
    plot(t, x, t, u(P)*x','LineWidth',1.0);
    makepretty('gcc_linear_bound2')

    figure(3)
    [t, x] = ode45(@(t, x) ...
	state(t, x, Gc(-1, -1)), [0 10], [1; 2]);
    plot(t, x, t, u(P)*x','LineWidth',1.0);
    makepretty('gcc_linear_bound3')
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
