function [P] = gcc_linear_bound()
    A = [-1, 0; 1, -2];
    B = [0; 1];
    Q = eye(length(A));
    R = 10;

    A_1 = @(r1) r1*[0, 0; 3.0, 0];
    A_2 = @(r2) r2*[0, 0; 0, 1.5];
    DA = @(r1, r2) A_1(r1) + A_2(r2);
    U = @(P, r1, r2, e) e*P + ...
	(1/e)*(DA(r1, r2)'*P*DA(r1, r2));

    e = 1;
    Pnot = care(A, B, Q, R);
    [P, fval] = fsolve(@(p) gare(p, A, B, Q, R, ...
	U(p, 1, 1, e)), Pnot);

    [P_t, fval_t] = fsolve(@(p_t) gcc_lyap(p_t, A, B, Q, R, ...
	U(p_t, 1, 1, e)), P);

    Gc = @(r1, r2) A + DA(r1, r2) - B*(R\(B'*P));
    makeplot(A, B, R, P, Pnot, Gc);
end

function f = gare(p, A, B, Q, R, U)
    P = reshape(p(1:length(A)^2), size(A));
    F = P*A + A'*P - P*B*(R\(B'*P)) + Q + U;
    f = F(:);
end

function f = gcc_lyap(p, A, B, Q, R, U)
    P = reshape(p(1:length(A)^2), size(A));
    K_t = -R\(B'*P);
    A_t = A + B*K_t;
    Q_t = Q + U;
    F = P*A_t + A_t'*P + Q_t + K_t'*R*K_t;
    f = F(:);
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
