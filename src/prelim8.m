A = [0, 1; 2, 1];
B = [0; 1];
C = [1, 1];
[K, s, e] = lqr(A, B, eye(2), 10, []);
step(ss(A - B*K, B, C, []));
grid on;
