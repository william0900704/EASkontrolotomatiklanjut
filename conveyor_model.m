function [A,B,C,D,plant] = conveyor_model()

A = [-29.63  -0.6085   1.181;
       2       0        0;
       0       1        0];
B = [4; 0; 0];
C = [0 0.0148 -5.2175];
D = 0;

plant = ss(A,B,C,D);
end
