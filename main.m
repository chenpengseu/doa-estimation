clear all; clc; close all;

%% system parameters
N = 16; % antenna number
d = 0.5; % normalized antenna distance
K = 6; % target number
sig_min_spacing = 100 / N; % deg, the minimal target spacing
spatial_angle_min = -60; % deg, the detection angle
spatial_angle_max = 60; % deg, the detection angle

% multiple trails to show the estimation performance
trail_num = 1e2;
RMSE_music = zeros(trail_num, 1);
RMSE_capon = zeros(trail_num, 1);
for idx_trail = 1:trail_num
    % genearte the random ground-truth DOA
    while (1) 
        if sig_min_spacing * K >= (spatial_angle_max - spatial_angle_min) * 0.7
            error('The targets are too close!');
        end 
        sig_angle = sort(rand(K, 1) * (spatial_angle_max - spatial_angle_min) + spatial_angle_min, 'ascend');
        if K > 1
            if (min(abs(sig_angle(2:end) - sig_angle(1: end - 1))) >= sig_min_spacing) 
                break;
            end 
        else
            break;
        end 
    end
    % signal and noise power
    SNR = -10; % dB
    Ps = ones(K, 1); % signal power 
    Pn = sum(Ps) / K / db2pow(SNR); 

    sig_len = 1e3; % signal length
    sig = bsxfun(@times, sqrt(Ps.' / 2), randn(sig_len, K) + 1j * randn(sig_len, K)); % target signals
    noise = sqrt(Pn / 2) * (randn(N, sig_len) + 1j * randn(N, sig_len));
    recv = get_steervec(N, d, deg2rad(sig_angle)) * sig.' + noise;

    ang_range = [max([spatial_angle_min - 10, -90]):0.01:min([90, spatial_angle_max + 10])].';
    ang_mat = get_steervec(N, d, deg2rad(ang_range));
    
    % music algorithm
    sp_music = music(recv, K, ang_mat);
    sp_music = sp_music / max(sp_music); 
    % get the estimated angle from the spectrum
    [est_ang_music, est_ang_index_music, RMSE_tmp] = get_estangle_from_spectrum(sp_music, ang_range, sig_angle, sig_min_spacing);
    RMSE_music(idx_trail) = RMSE_tmp;

    % capon algorithm
    sp_capon = capon(recv, K, ang_mat);
    sp_capon = sp_capon / max(sp_capon);
    % get the estimated angle from the spectrum
    [est_ang_capon, est_ang_index_capon, RMSE_tmp] = get_estangle_from_spectrum(sp_capon, ang_range, sig_angle, sig_min_spacing);
    RMSE_capon(idx_trail) = RMSE_tmp;

    if idx_trail==1
        % show the spectrum
        sp_music = pow2db(sp_music);
        sp_capon = pow2db(sp_capon);
        figure; plot(ang_range, sp_music);
        hold on;
        stem(est_ang_music, sp_music(est_ang_index_music), 'x', 'LineStyle', 'none');
        stem(sig_angle, pow2db(Ps / max(Ps)), '-o', 'BaseValue', -80);
        plot(ang_range, sp_capon); 
        stem(est_ang_capon, sp_capon(est_ang_index_capon), 'x', 'LineStyle', 'none');

        legend('MUSIC spectrum', 'Estimated angles (MUSIC)', 'Ground-truth angles',...
             'Capon spectrum', 'Estimated angles (Capon)')
        grid on;
        drawnow;
    end
end
fprintf('RMSE (MUSIC): %.4g deg\n', sqrt(sum(abs(RMSE_music).^2) / length(RMSE_music)));
fprintf('RMSE (Capon): %.4g deg\n', sqrt(sum(abs(RMSE_capon).^2) / length(RMSE_capon)));
