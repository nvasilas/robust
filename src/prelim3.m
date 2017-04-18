sys1 = tf([4], [1, 2, 3]);
sys2 = 1;
sys = feedback(sys1, sys2);
