function [P] = gcc_linear_bound()
    A = [0, 1; -2, 1];
    B = [0; 1];
    Q = eye(length(A));
    R = 10;

    A_1 = @(r1) r1*[0, 0; 0.5, 0];
    A_2 = @(r2) r2*[0, 0; 0, 0.8];
    DA = @(r1, r2) A_1(r1) + A_2(r2);
    U = @(P, r1, r2, e) e*P + ...
	(1/e)*(DA(r1, r2)'*P*DA(r1, r2));

    e = 1;
    Pnot = care(A, B, Q, R);
    [P, fval] = fsolve(@(p) gare(p, A, B, Q, R, ...
	U(p, 1, 1, e)), Pnot);

    [P_t, fval_t] = fsolve(@(p_t) gcc_lyap(p_t, A, B, Q, R, ...
	U(p_t, 1, 1, e)), P);

    Gc = @(r1, r2, p1) A + DA(r1, r2) - B*(R\(B'*P));
    makeplot(A, B, R, Pnot, Gc);
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

function makeplot(A, B, R, Pnot, Gc)
    close all;
    figure(1)
    [t, x] = ode45(@(t, x) state(t, x, A - B*(R\(B'*Pnot))), [0 15], [1; 2]);
    plot(t, x)
    title('Initial System')
    grid on;

    figure(2)
    [t, x] = ode45(@(t, x) state(t, x, Gc(1, 1, 1)), [0 15], [1; 2]);
    plot(t, x)
    title('Uncertain System r1 = 1, r2 = 1, p1 = 1')
    grid on;

    figure(3)
    [t, x] = ode45(@(t, x) state(t, x, Gc(-1, -1, -0.5)), [0 20], [1; 2]);
    plot(t, x)
    title('Uncertain System r1 = -1, r2 = -1, p1 = -0.5')
    grid on;
end
