function TEX_pp(p, K, Gc, state);
    addpath('~/downloads/matlab2tikz-master/src/');

    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, 0, 0)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    makepretty('./pretty_fig/pp_ex1')

    figure(2)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, 1, 1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    makepretty('./pretty_fig/pp_ex2')

    figure(3)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, -1, -1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    makepretty('./pretty_fig/pp_ex3')
end

function makepretty(string)
    xlabel('Time (seconds)','interpreter','latex');
    ylabel('Amplitude','interpreter','latex');
    grid on;
    l = legend(['$x_1(t)$'], ['$x_2(t)$'], ...
        ['$x_3(t)$'], ['$u(t)$']);
    set(l,'Interpreter','Latex');
    cleanfigure;
    matlab2tikz(strcat(string, '.tex'));
end
