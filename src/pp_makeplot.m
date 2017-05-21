function pp_makeplot(p)
    DA = @(r1, r2) ...
        r1*p.D(:,1)*p.E(1,:) ...
        + r2*p.D(:,2)*p.E(2,:);
    Gc = @(r1, r2) ...
        p.A + DA(r1, r2) - p.B*(p.R\(p.B'*p.P_f));
    K = -p.R\(p.B'*p.P_f);

    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(1, 1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    title('Uncertain System r1 = 1, r2 = 1')
    legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'u(t)');
    grid on;

    figure(2)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(-1, -1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    title('Uncertain System r1 = -1, r2 = -1')
    legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'u(t)');
    grid on;

    maketable(p, Gc);
end

function dxdt = state(~, x, Gc)
    dxdt = Gc*x;
end

function maketable(p, Gc)
    fprintf('r1 = -1, r2 = -1 closed loop eigenvalues\n');
    eig(Gc(-1, -1))

    fprintf('r1 = 1, r2 = 1 closed loop eigenvalues\n');
    eig(Gc(1, 1))

    fprintf('r1 = 1, r2 = -1 closed loop eigenvalues\n');
    eig(Gc(1, -1))

    fprintf('r1 = 0, r2 = -1 closed loop eigenvalues\n');
    eig(Gc(0, -1))

    fprintf('r1 = 1, r2 = 0 closed loop eigenvalues\n');
    eig(Gc(1, 0))
end
