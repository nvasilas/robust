function pp_makeplot(p, string)
    K = -p.R\(p.B'*p.P_f);

    if strcmp(string, 'TEX')
        addpath('./pretty_fig');
        TEX_pp(p, K, @Gc, @state);
        return
    end

    close all;
    figure(1)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, 0, 0)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    title('Nominal System r1 = 0, r2 = 0')
    legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'u(t)');
    grid on;

    figure(2)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, 1, 1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    title('Uncertain System r1 = 1, r2 = 1')
    legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'u(t)');
    grid on;

    figure(3)
    [t, x] = ode45(@(t, x) ...
        state(t, x, Gc(p, -1, -1)), [0 6], [3; -2; 2]);
    plot(t, x, t, K*x','LineWidth',1.0);
    title('Uncertain System r1 = -1, r2 = -1')
    legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'u(t)');
    grid on;

    maketable(p);
end

function closed = Gc(p, r1, r2)
    if (r1 > -eps) && (r1 < eps) ...
            && (r2 > -eps) && (r2 < eps)
        [closed, ~, error_flag] = pp_nominal(p);
        if error_flag
            fprintf('Error H, no distinct eigenvalues\n');
            return
        end
        return
    end
    DA = @(r1, r2) ...
        r1*p.D(:,1)*p.E(1,:) ...
        + r2*p.D(:,2)*p.E(2,:);
    closed = p.A + DA(r1, r2) - p.B*(p.R\(p.B'*p.P_f));
end

function dxdt = state(~, x, K)
    dxdt = K*x;
end

function print_msg(lambda)
    fprintf('lambda_1 = %s\n', num2str(lambda(1)));
    fprintf('lambda_2 = %s\n', num2str(lambda(2)));
    fprintf('lambda_3 = %s\n', num2str(lambda(3)));
end

function maketable(p)
    fprintf('closed loop eigenvalues for r1 = 0, r2 = 0\n');
    print_msg(eig(Gc(p, 0, 0)));

    fprintf('\n closed loop eigenvalues for r1 = -1, r2 = -1\n');
    print_msg(eig(Gc(p, -1, -1)));

    fprintf('\n closed loop eigenvalues for r1 = 1, r2 = 1\n');
    print_msg(eig(Gc(p, 1, 1)));

    fprintf('\n closed loop eigenvalues for r1 = 1, r2 = -1\n');
    print_msg(eig(Gc(p, 1, -1)));

    fprintf('\n closed loop eigenvalues for r1 = 0, r2 = -1\n');
    print_msg(eig(Gc(p, 0, -1)));

    fprintf('\n closed loop eigenvalues for r1 = 1, r2 = 0\n');
    print_msg(eig(Gc(p, 1, 0)));
end
