sys = tf([1], [1, 1.54, 1]);
[u, t] = gensig('square', 20, 50, 0.1);
lsim(sys, u, t);
grid on;
