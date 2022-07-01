clear all; clc; close all;

%% system parameters
N = 16; % antenna number
d = 0.5; % normalized antenna distance
K = 3; % target number
sig_angle = [-30, 0, 20].'; % target angle
sig_len = 1e4; % signal length
Ps = ones(K, 1); % signal power
sig = bsxfun(@times, sqrt(Ps.' / 2), randn(sig_len, K) + 1j * randn(sig_len, K)); % target signals
SNR = 20; % dB
Pn = sum(Ps) / K / db2pow(SNR);
noise = sqrt(Pn / 2) * (randn(N, sig_len) + 1j * randn(N, sig_len));
recv = get_steervec(N, d, deg2rad(sig_angle)) * sig.' + noise;

ang_range = [-70:0.01:70].';
ang_mat = get_steervec(N, d, deg2rad(ang_range));
sp_music = music(recv, K, ang_mat);

sp_music = pow2db(sp_music / max(sp_music));
figure; plot(ang_range, sp_music);
hold on;
stem(sig_angle, pow2db(Ps / max(Ps)), '-o', 'BaseValue', -80);
legend('MUSIC', 'Ground-truth angles')
grid on;

