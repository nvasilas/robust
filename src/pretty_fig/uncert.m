%% Define the matrix A
A = [0 1 0; 0 0 1; -1 -3 2];
eig(A);
%Since the system is in cannonical controllability form B=[0 0 1]'
B = [0 0 1]';
C = [1 1 1];
D = 0;
%% Apply LQR with Q=eye(3) and R=1
Q = eye(3);
R = 1;

K = lqr(A, B, Q, R);
%The state matrix for the stabilized (Ac) system is A-BK
Ac = A - B*K;
eig(Ac);

%% Define the uncertainties
%Suppose uncertainties affect the elements in the positions
%(3,2) and (3,3) of the initial matrix A
%Thus the matrices that state the positions of each
%uncertainties in A are

A1 = [0 0 0;0 0 0; 0 1 0];
A2 = [0 0 0;0 0 0; 0 0 1];

%Also define the range of each uncertainty, namely r1 and r2.
r1 = 1;
r2 = 1;
%NOTE:
%Here the range does not affect the matrices A1 and A2, but
%this will now always be the case

%% Stability through uncertainties check
%Now check if K that was calculated for the matrix A,
%stabilizes all the uncertain systems that arise

figure(1)
cleanfigure;
hold on;
for i = -10:1:10
    for j = -10:1:10
        A_uncertain = A + A1*i + A2*j - B*K;
        eigvals = eig(A_uncertain)';
        %using max of eigenvalue vector > 0 means that at
	%least one of the eigenvalues are positive and
	%thus the system is unstable
        if max(real(eigvals)) >= 0
            plot(i, j, '*r')
        else
            plot(i, j, '+g')
        end
    end
end
grid on;
xlabel('Uncertainty $r_1$','interpreter','latex');
ylabel('Uncertainty $r_2$','interpreter','latex');
matlab2tikz('uncert_r.tex');

%% Plot the step responces for one of the stable systems
%that arise from certain values of uncertainties.
% Also plot the initial stabilized system and
%compare it with the uncertain one.
figure(2)
hold on;
cleanfigure;
sys1 = ss(Ac, B, C, D);
sys2 = ss(A + A1*2 + A2*2 - B*K, B, C, D);
[y1, ~] = step(sys1);
[y2, t]  = step(sys2);
plot(t, y1, t, y2, 'LineWidth', 1.0);
plot([t(1) t(end)], [0.707106 0.707106], 'k:', 'LineWidth', 1.0);
xlim([0 t(end)]);
grid on;
l = legend('Initial Stabilized System','Uncertain Stable System');
set(l,'Interpreter','Latex');
xlabel('Time (seconds)','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
%matlab2tikz('uncert_step.tex');
