clear; clc; close all;

[A,B,C,D,plant] = conveyor_model();
n = size(A,1);

% Bagian LQR
Q_lqr = diag([1 1 10000]);
R_lqr = 0.1;

[K_lqr,~,~] = lqr(A,B,Q_lqr,R_lqr);

Acl_lqr = A - B*K_lqr;

% Feedforward gain untuk tracking unit step
Nbar = -1/(C*(Acl_lqr\B));

% Bagian Kalman Filter

G = eye(n);

Qn = diag([1e-4 1e-4 1e-4]);   % process-noise covariance
Rn = 1e-3;                      % measurement-noise covariance
Nn = zeros(n,1);

[L_kf,Pkf,E_kf] = lqe(A,G,C,Qn,Rn,Nn);

% Sistem Closed-Loop Augmented LQG

A_lqg = [A,            -B*K_lqr;
         L_kf*C, A - B*K_lqr - L_kf*C];

B_lqg = [B*Nbar;
         B*Nbar];

C_lqg = [C zeros(1,n)];
D_lqg = 0;

sys_lqg = ss(A_lqg, B_lqg, C_lqg, D_lqg);

% Analisis Kestabilan Pole
pole_open_loop = eig(A);
pole_regulator = eig(Acl_lqr);
pole_estimator = eig(A - L_kf*C);
pole_lqg = eig(A_lqg);

is_open_loop_stable = all(real(pole_open_loop) < 0);
is_regulator_stable = all(real(pole_regulator) < 0);
is_estimator_stable = all(real(pole_estimator) < 0);
is_lqg_stable = all(real(pole_lqg) < 0);

fprintf('ANALISIS KESTABILAN POLE \n');

fprintf('\nPole open-loop:\n');
disp(pole_open_loop);

if is_open_loop_stable
    fprintf('Sistem open-loop stabil karena semua pole berada di kiri bidang-s.\n');
else
    fprintf('Sistem open-loop tidak stabil karena terdapat pole dengan bagian real positif.\n');
end

fprintf('\nPole regulator LQR A-BK:\n');
disp(pole_regulator);

if is_regulator_stable
    fprintf('Regulator LQR stabil karena semua pole A-BK berada di kiri bidang-s.\n');
else
    fprintf('Regulator LQR tidak stabil.\n');
end

fprintf('\nPole estimator Kalman A-LC:\n');
disp(pole_estimator);

if is_estimator_stable
    fprintf('Estimator Kalman stabil karena semua pole A-LC berada di kiri bidang-s.\n');
else
    fprintf('Estimator Kalman tidak stabil.\n');
end

fprintf('\nPole closed-loop augmented LQG:\n');
disp(pole_lqg);

if is_lqg_stable
    fprintf('Sistem closed-loop LQG stabil karena semua pole augmented berada di kiri bidang-s.\n\n');
else
    fprintf('Sistem closed-loop LQG tidak stabil karena masih terdapat pole dengan bagian real positif.\n\n');
end

% Analisis Kestabilan Lyapunov

Q_lyap = eye(size(A_lqg));

P_lyap = lyap(A_lqg', Q_lyap);

eig_P = eig(P_lyap);
lyap_residual = A_lqg'*P_lyap + P_lyap*A_lqg;

is_P_positive_definite = all(eig_P > 0);

fprintf('=== ANALISIS KESTABILAN LYAPUNOV LQG ===\n');

fprintf('Matriks P hasil persamaan Lyapunov:\n');
disp(P_lyap);

fprintf('Eigenvalue matriks P:\n');
disp(eig_P);

fprintf('Residual A_lqg''P + P*A_lqg:\n');
disp(lyap_residual);

if is_P_positive_definite && is_lqg_stable
    fprintf('Sistem closed-loop LQG stabil secara Lyapunov karena P positive definite.\n\n');
else
    fprintf('Sistem closed-loop LQG tidak memenuhi syarat kestabilan Lyapunov.\n\n');
end

% Simulasi Respons Step Tanpa Noise
Tend = 3;
t = 0:0.001:Tend;
r = ones(size(t));

[y_lqg,t_lqg,z_lqg] = lsim(sys_lqg, r, t);

xhat = z_lqg(:,n+1:end);

u_lqg = (-K_lqr*xhat.' + Nbar*r).';

info_lqg_no_noise = stepinfo(y_lqg, t_lqg, 1);

% Simulasi LQG dengan Noise
rng(10);

w = 0.01*randn(length(t),n);      % process noise
v = 0.01*randn(length(t),1);      % measurement noise

% Input sistem augmented: [r, w, v]
B_aug = [B*Nbar, G,          zeros(n,1);
         B*Nbar, zeros(n,n), L_kf];

C_aug = [C zeros(1,n)];
D_aug = zeros(1,1+n+1);

sys_lqg_noise = ss(A_lqg, B_aug, C_aug, D_aug);

U_aug = [r(:), w, v];

[y_lqg_noise,~,z_noise] = lsim(sys_lqg_noise, U_aug, t);

% Output sensor yang noisy
y_meas_noise = y_lqg_noise + v;

% State estimasi saat noise
xhat_noise = z_noise(:,n+1:end);

u_lqg_noise = (-K_lqr*xhat_noise.' + Nbar*r).';

% Analisis Performansi LQG Tanpa dan Dengan Noise
% ===============================

% Performansi output aktual saat noise
info_lqg_noise_true = stepinfo(y_lqg_noise, t, 1);

% Performansi output terukur yang terkena noise
info_lqg_noise_meas = stepinfo(y_meas_noise, t, 1);

% Error tracking terhadap referensi
e_lqg_no_noise = r(:) - y_lqg;
e_lqg_noise_true = r(:) - y_lqg_noise;
e_lqg_noise_meas = r(:) - y_meas_noise;

% RMSE tracking
RMSE_no_noise = sqrt(mean(e_lqg_no_noise.^2));
RMSE_noise_true = sqrt(mean(e_lqg_noise_true.^2));
RMSE_noise_meas = sqrt(mean(e_lqg_noise_meas.^2));

% Tabel performansi
Performance_LQG = table( ...
    [info_lqg_no_noise.RiseTime; info_lqg_noise_true.RiseTime; info_lqg_noise_meas.RiseTime], ...
    [info_lqg_no_noise.SettlingTime; info_lqg_noise_true.SettlingTime; info_lqg_noise_meas.SettlingTime], ...
    [info_lqg_no_noise.Overshoot; info_lqg_noise_true.Overshoot; info_lqg_noise_meas.Overshoot], ...
    [info_lqg_no_noise.PeakTime; info_lqg_noise_true.PeakTime; info_lqg_noise_meas.PeakTime], ...
    [RMSE_no_noise; RMSE_noise_true; RMSE_noise_meas], ...
    'VariableNames', {'RiseTime_s','SettlingTime_s','Overshoot_percent','PeakTime_s','RMSE'}, ...
    'RowNames', {'LQG_Tanpa_Noise','LQG_Dengan_Noise_Output_Aktual','LQG_Dengan_Noise_Output_Terukur'} );

% Output Hasil Rancangan
fprintf('RANCANGAN LQG \n');

fprintf('Q_lqr =\n');
disp(Q_lqr);

fprintf('R_lqr = %.4g\n', R_lqr);

fprintf('K_lqr =\n');
disp(K_lqr);

fprintf('Nbar = %.6g\n', Nbar);

fprintf('Qn =\n');
disp(Qn);

fprintf('Rn = %.6g\n', Rn);

fprintf('Kalman gain L =\n');
disp(L_kf);

fprintf('Step info LQG tanpa noise:\n');
disp(info_lqg_no_noise);

fprintf('Step info LQG dengan noise, output aktual:\n');
disp(info_lqg_noise_true);

fprintf('Step info LQG dengan noise, output terukur:\n');
disp(info_lqg_noise_meas);

fprintf('=== TABEL PERFORMANSI LQG ===\n');
disp(Performance_LQG);

writetable(Performance_LQG, 'performance_lqg_noise.csv', 'WriteRowNames', true);

%% ==============================
% Plot Respons Step LQG Tanpa Noise
% ===============================
figure('Color','w');
plot(t_lqg, y_lqg, 'LineWidth', 1.8); hold on;
yline(1,'--','Reference');

grid on;
xlabel('Time (s)');
ylabel('Belt speed response (normalized)');
title('Respons Step Closed-Loop LQG Tanpa Noise');

% Plot LQG dengan Noise Pengukuran
figure('Color','w');
plot(t, y_meas_noise, ':', 'LineWidth', 1.0); hold on;
plot(t, y_lqg_noise, 'LineWidth', 1.8);
yline(1,'--','Reference');

grid on;
xlabel('Time (s)');
ylabel('Belt speed response (normalized)');
title('Respons LQG dengan Noise Pengukuran');
legend('Measured output + noise','Output aktual sistem','Reference','Location','best');

% Plot Perbandingan Tanpa Noise dan Dengan Noise
figure('Color','w');
plot(t_lqg, y_lqg, 'LineWidth', 1.8); hold on;
plot(t, y_lqg_noise, '--', 'LineWidth', 1.8);
plot(t, y_meas_noise, ':', 'LineWidth', 1.0);
yline(1,'--','Reference');

grid on;
xlabel('Time (s)');
ylabel('Belt speed response (normalized)');
title('Perbandingan Respons LQG Tanpa dan Dengan Noise');
legend('LQG tanpa noise', ...
       'LQG dengan noise, output aktual', ...
       'LQG dengan noise, output terukur', ...
       'Reference', ...
       'Location','best');

% Plot Sinyal Kontrol
figure('Color','w');
plot(t_lqg, u_lqg, 'LineWidth', 1.8); hold on;
plot(t, u_lqg_noise, '--', 'LineWidth', 1.5);

grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Sinyal Kontrol LQG Tanpa dan Dengan Noise');
legend('LQG tanpa noise','LQG dengan noise','Location','best');