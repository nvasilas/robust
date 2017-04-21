sys = tf([1], [1, 1.54, 1]);
[u, t] = gensig('square', 20, 50, 0.1);
y = lsim(sys, u, t);
hold on;
plot(t(1:end-1), y(1:end-1), 'LineWidth', 1.0);
plot(t(1:end-1), u(1:end-1),'k', 'LineWidth', 1.0);
xlabel('Time (seconds)','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
grid on;
matlab2tikz('prelim7.tex');