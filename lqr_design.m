clear; clc; close all;

[A,B,C,D,plant] = conveyor_model();

% Bobot LQR
Q = diag([1 1 10000]);
R = 0.1;

[K_lqr,S_lqr,e_lqr] = lqr(A,B,Q,R);

Acl = A - B*K_lqr;

% Feedforward gain untuk tracking unit step
Nbar_lqr = -1/(C*(Acl\B));

sys_lqr = ss(Acl, B*Nbar_lqr, C, D);

% Analisis Kestabilan Pole
pole_open_loop = eig(A);
pole_closed_loop = eig(Acl);

is_open_loop_stable = all(real(pole_open_loop) < 0);
is_closed_loop_stable = all(real(pole_closed_loop) < 0);

fprintf('ANALISIS KESTABILAN POLE \n');

fprintf('Pole open-loop:\n');
disp(pole_open_loop);

if is_open_loop_stable
    fprintf('Sistem open-loop stabil karena semua pole berada di kiri bidang-s.\n\n');
else
    fprintf('Sistem open-loop tidak stabil karena terdapat pole dengan bagian real positif.\n\n');
end

fprintf('Pole closed-loop LQR:\n');
disp(pole_closed_loop);

if is_closed_loop_stable
    fprintf('Sistem closed-loop LQR stabil karena semua pole berada di kiri bidang-s.\n\n');
else
    fprintf('Sistem closed-loop LQR tidak stabil karena masih terdapat pole dengan bagian real positif.\n\n');
end

% Analisis Kestabilan Lyapunov

Q_lyap = eye(size(Acl));

P_lyap = lyap(Acl', Q_lyap);

eig_P = eig(P_lyap);
lyap_residual = Acl'*P_lyap + P_lyap*Acl;

is_P_positive_definite = all(eig_P > 0);

fprintf('ANALISIS KESTABILAN LYAPUNOV \n');

fprintf('Matriks P hasil persamaan Lyapunov:\n');
disp(P_lyap);

fprintf('Eigenvalue matriks P:\n');
disp(eig_P);

fprintf('Residual Acl''P + P*Acl:\n');
disp(lyap_residual);

if is_P_positive_definite && is_closed_loop_stable
    fprintf('Sistem closed-loop LQR stabil secara Lyapunov karena P positive definite.\n\n');
else
    fprintf('Sistem closed-loop LQR tidak stabil dengan syarat Lyapunov.\n\n');
end

% Simulasi Respon Step
info_lqr = stepinfo(sys_lqr);

t = 0:0.001:3;
r = ones(size(t));

[y_lqr,t_lqr,x_lqr] = lsim(sys_lqr, r, t);

u_lqr = -K_lqr*x_lqr.' + Nbar_lqr*r;
u_lqr = u_lqr.';

fprintf(' RANCANGAN LQR \n');
fprintf('Q =\n'); disp(Q);
fprintf('R = %.4g\n', R);
fprintf('K_lqr =\n'); disp(K_lqr);
fprintf('Nbar_lqr = %.6g\n', Nbar_lqr);
fprintf('Step info LQR:\n'); disp(info_lqr);

% Plot Respon Step
figure('Color','w');
plot(t_lqr, y_lqr, 'LineWidth', 1.8); hold on;
yline(1,'--','Reference');

grid on;
xlabel('Time (s)');
ylabel('Belt speed response (normalized)');
title('Respons Step Closed-Loop LQR');

% Plot Sinyal Kontrol
figure('Color','w');
plot(t_lqr, u_lqr, 'LineWidth', 1.8);

grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Sinyal Kontrol LQR');