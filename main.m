% This is the main function to test the estimation methods and show the DOA estimation performance.

clear all; clc; close all;

%% system parameters
N = 16; % antenna number
d = 0.5; % normalized antenna distance
K = 3; % target number
sig_min_spacing = 100 / N*2; % deg, the minimal target spacing
spatial_angle_min = -60; % deg, the detection angle
spatial_angle_max = 60; % deg, the detection angle
SNR = 10; % dB
ang_grid = 0.01;
dic_grid = 0.1;
is_offgrid = 1;

% dictionary matrix
dic_range = [max([spatial_angle_min - 10, -90]):dic_grid:min([90, spatial_angle_max + 10])].';
dic_mat = get_steervec(N, d, deg2rad(dic_range));
ang_range = [max([spatial_angle_min - 10, -90]):ang_grid:min([90, spatial_angle_max + 10])].';
ang_mat = get_steervec(N, d, deg2rad(ang_range));

% multiple trails to show the estimation performance
trail_num = 10;
RMSE_music = zeros(trail_num, 1);
RMSE_capon = zeros(trail_num, 1);
RMSE_somp = zeros(trail_num, 1);
t_capon = 0;
t_music = 0;
t_somp = 0;
for idx_trail = 1:trail_num
    % genearte the random ground-truth DOA
    while (1) 
        if sig_min_spacing * K >= (spatial_angle_max - spatial_angle_min) * 0.7
            error('The targets are too close!');
        end 


        if is_offgrid
            sig_angle = sort(rand(K, 1) * (spatial_angle_max - spatial_angle_min) + spatial_angle_min, 'ascend');
        else
            sig_idx = randi([min(find(spatial_angle_min<dic_range)), max(find(spatial_angle_max>dic_range))], K, 1);
            sig_angle = sort(dic_range(sig_idx), 'ascend');
        end
        if K > 1
            if (min(abs(sig_angle(2:end) - sig_angle(1: end - 1))) >= sig_min_spacing) 
                break;
            end 
        else
            break;
        end 
    end
    % signal and noise power
    Ps = ones(K, 1); % signal power 
    Pn = sum(Ps) / K / db2pow(SNR); 

    % generate the signals
    sig_len = 1e3; % signal length
    sig = bsxfun(@times, sqrt(Ps.' / 2), randn(sig_len, K) + 1j * randn(sig_len, K)); % target signals
    noise = sqrt(Pn / 2) * (randn(N, sig_len) + 1j * randn(N, sig_len));
    recv = get_steervec(N, d, deg2rad(sig_angle)) * sig.' + noise;

    
    %% music algorithm
    tic;
    sp_music = music(recv, K, ang_mat);
    t_music = t_music+toc;
    sp_music = sp_music / max(sp_music); 
    % get the estimated angle from the spectrum
    [est_ang_music, est_ang_index_music, RMSE_tmp] = get_estangle_from_spectrum(sp_music, ang_range, sig_angle, sig_min_spacing);
    RMSE_music(idx_trail) = RMSE_tmp;

    %% capon algorithm
    tic
    sp_capon = capon(recv, K, ang_mat);
    t_capon = t_capon+toc;
    sp_capon = sp_capon / max(sp_capon);
    % get the estimated angle from the spectrum
    [est_ang_capon, est_ang_index_capon, RMSE_tmp] = get_estangle_from_spectrum(sp_capon, ang_range, sig_angle, sig_min_spacing);
    RMSE_capon(idx_trail) = RMSE_tmp;

    %% SOMP
    tic;
    sp_somp = somp(recv, K, dic_mat);
    t_somp = t_somp+toc;
    sp_somp = sp_somp / max(sp_somp);  
    % get the estimated angle from the spectrum
    [est_ang_somp, est_ang_index_somp, RMSE_tmp] = get_estangle_from_spectrum(sp_somp, dic_range, sig_angle, sig_min_spacing);
    RMSE_somp(idx_trail) = RMSE_tmp;

    if idx_trail==1
        last_time = (t_somp+t_music+t_capon)*trail_num;
        cur_time = 0;
        fprintf('Waiting for %.2g s...\n', last_time);
    end
    if idx_trail==1
        % show the spectrum
        sp_music = pow2db(sp_music);
        sp_capon = pow2db(sp_capon);
        sp_somp = pow2db(sp_somp);
        figure; 
        leg_str = {};
        stem(sig_angle, pow2db(Ps / max(Ps)), '-o', 'BaseValue', -80, 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{1} = 'Ground-truth angles';
        hold on;
        plot(ang_range, sp_music, 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{2} = 'MUSIC spectrum'; 
        stem(est_ang_music, sp_music(est_ang_index_music), 'x', 'LineStyle', 'none', 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{3} = '';
        plot(ang_range, sp_capon, 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{4} = 'Capon spectrum';
        stem(est_ang_capon, sp_capon(est_ang_index_capon), 'x', 'LineStyle', 'none', 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{5} = ''; 
        stem(dic_range, sp_somp, '-s', 'BaseValue', -80, 'LineWidth', 1, 'MarkerSize', 8);
        leg_str{6} = 'SOMP'; 
        legend(leg_str);
        grid on;
        drawnow;
    end

end
fprintf('\n**** DOA Estimation Methods ***** \n');
fprintf('Parameters:\n');
fprintf('Antenna number: %d\n', N);
fprintf('Antenna distance: %.2g wavelength\n', d);
fprintf('Target number: %d\n', K);
fprintf('Minimual target spacing: %.2g deg\n', sig_min_spacing); 
fprintf('Grid size (non-grid method): %.2g deg\n', ang_grid);
fprintf('Grid size (grid method): %.2g deg\n', dic_grid);
fprintf('SNR: %.2g dB\n', SNR);
fprintf('Targets off-grid: %d \n', is_offgrid);
fprintf('-----------\n')
fprintf('Performance:\n');
fprintf('RMSE (MUSIC): %.4g deg, Time: %.4g s\n', sqrt(sum(abs(RMSE_music).^2) / length(RMSE_music)), t_music/trail_num);
fprintf('RMSE (Capon): %.4g deg, Time: %.4g s\n', sqrt(sum(abs(RMSE_capon).^2) / length(RMSE_capon)), t_capon / trail_num);
fprintf('RMSE (SOMP): %.4g deg, Time: %.4g s\n', sqrt(sum(abs(RMSE_somp).^2) / length(RMSE_somp)), t_somp / trail_num);
fprintf('*************************************\n');

