% This is the main function to test the estimation methods 
% and show the DOA estimation performance.
function get_doa_snr(method)
    method_idx = 0;
    switch method
    case 'music'
        method_idx = 1;
    case 'capon'
        method_idx = 2;
    case 'somp'
        method_idx = 3;
    case 'esprit'
        method_idx = 4;
    case 'root-music'
        method_idx = 5;
    otherwise
        error('Method Error!');
    end

    [param] = init_system();
 
    % dictionary matrix for somp method
    dic_grid = 0.1; 
    dic_range = [max([param.doa_min - 10, -90]):dic_grid:min([90, param.doa_max + 10])].';
    dic_mat = get_steervec(param.N, param.d, deg2rad(dic_range));

    ang_grid = 0.01;
    doa_range = [max([param.doa_min - 10, -90]):ang_grid:min([90, param.doa_max + 10])].';
    ang_mat = get_steervec(param.N, param.d, deg2rad(doa_range));

    % multiple trails to show the estimation performance
    sig_len = 1e3; % signal length
    trail_num = 1e3;
    SNR_range = [-20:10:40].'-10*log10(sig_len);
    t = 0;
    RMSE = zeros(trail_num, length(SNR_range));
    for idx_SNR = 1:length(SNR_range)
        SNR = SNR_range(idx_SNR); % dB
        for idx_trail = 1:trail_num
            % genearte the random ground-truth DOA
            doa = generate_doa(param); 

            recv = generate_recv(param, sig_len, doa, SNR);
            switch method
            case 'music'
                %% music algorithm
                tic; 
                sp_music = music(recv, param.K, ang_mat);
                t = t+toc;
                sp_music = sp_music / max(sp_music); 
                % get the estimated angle from the spectrum
                [est_ang_music, est_ang_index_music, RMSE_tmp] = get_doa_from_spectrum(sp_music, doa_range, doa, param.doa_min_spacing);
                RMSE(idx_trail, idx_SNR) = RMSE_tmp;
            case 'capon'
                %% capon algorithm
                tic
                sp_capon = capon(recv, param.K, ang_mat);
                t = t+toc;
                sp_capon = sp_capon / max(sp_capon);
                % get the estimated angle from the spectrum
                [est_ang_capon, est_ang_index_capon, RMSE_tmp] = get_doa_from_spectrum(sp_capon, doa_range, doa, param.doa_min_spacing);
                RMSE(idx_trail, idx_SNR) = RMSE_tmp;
            case 'somp'
                %% SOMP
                tic;
                sp_somp = somp(recv, param.K, dic_mat);
                t = t+toc;
                sp_somp = sp_somp / max(sp_somp);  
                % get the estimated angle from the spectrum
                [est_ang_somp, est_ang_index_somp, RMSE_tmp] = get_doa_from_spectrum(sp_somp, dic_range, doa, param.doa_min_spacing);
                RMSE(idx_trail, idx_SNR) = RMSE_tmp;
            case 'esprit'
                %% ESPRIT 
                tic; 
                est_doa = ESPRIT(recv, param.K, param.d);
                t = t + toc; 
                rmse = get_rmse(doa, est_doa);
                RMSE(idx_trail, idx_SNR) = rmse; 
            case 'root-music'
                %% Root-MUSIC algorithm
                tic;
                est_doa = root_music(recv, param.K, param.d);
                t = t + toc;
                rmse = get_rmse(doa, est_doa);
                RMSE(idx_trail, idx_SNR) = rmse;
            otherwise
                error('Method Error!'); 
            end 
    end
end
RMSE_SNR= vec(sqrt(sum(abs(RMSE).^2) / size(RMSE, 1)));

% save the files
if ~exist('SNR.mat')
    save('SNR.mat', 'SNR_range');
    SNR_tmp = SNR_range;
else
    SNR_file = load('SNR.mat');
    SNR_tmp = SNR_file.SNR_range;
end
if norm(SNR_tmp-SNR_range)>1e-6 || ~exist('RMSE_SNR.mat')
    RMSE_SNR_all = [];
else
    RMSE_SNR_file = load('RMSE_SNR.mat');
    RMSE_SNR_all = RMSE_SNR_file.RMSE_SNR_all;
    if size(RMSE_SNR_all, 1)~=length(RMSE_SNR)
        RMSE_SNR_all = [];
    end
end

RMSE_SNR_all(:, method_idx) = RMSE_SNR;
save('RMSE_SNR.mat', 'RMSE_SNR_all');

if ~exist('t.mat')
    t_all = [];
else
    t_file = load('t.mat');
    t_all = t_file.t_all;
end
t_all(method_idx) = t/length(SNR_range)/trail_num;
save('t.mat', 't_all');

% show results
fprintf('\n******* %s *******\n', method)
fprintf('Time(s): %.4g\n', t_all(method_idx));
for idx_SNR = 1:length(SNR_range)
    SNR = SNR_range(idx_SNR); % dB
    fprintf('SNR(dB): %.2f,\t RMSE(deg): %.4g\n', SNR, RMSE_SNR(idx_SNR));
end

legend_str = {'MUSIC', 'Capon', 'SOMP', 'ESPRIT', 'Root\_MUSIC'};
h=figure; 
semilogy(SNR_range, RMSE_SNR_all, 'LineWidth', 2)
legend(legend_str{1:size(RMSE_SNR_all,2)});
grid on;
set(get(gca, 'XLabel'), 'String', 'SNR(dB)');
set(get(gca, 'YLabel'), 'String', 'RMSE (deg)');
drawnow;
savefig(h, './figures/SNR-RMSE');
 
