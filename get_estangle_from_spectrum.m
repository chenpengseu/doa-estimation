function [est_ang, est_ang_index, RMSE] = get_estangle_from_spectrum(sp, ang_range, sig_angle, sig_min_spacing)
    % get the estimated angle from the spectrum
    [pks, locs] = findpeaks(sp);
    [pks, sort_idx] = sort(pks, 'descend');
    locs = locs(sort_idx);
    ang_pks = ang_range(locs);
    valid_pks = 0;
    valid_pkx_idx = 1;
    if length(ang_pks) > 0
        valid_pks = [];
        valid_pkx_idx = [];
        for idx = 1:length(ang_pks)

            if idx == 1
                valid_pks(1) = ang_pks(idx);
                valid_pkx_idx(1) = locs(idx);
            else

                if min(abs(ang_pks(idx) -valid_pks)) > sig_min_spacing
                    valid_pks(end + 1) = ang_pks(idx);
                    valid_pkx_idx(end + 1) = locs(idx);
                end

            end

        end

    end

    [~, min_idx] = min(abs(vec(sig_angle) - vecH(valid_pks)).');
    est_ang = vec(valid_pks(min_idx));
    est_ang_index = vec(valid_pkx_idx(min_idx));
    RMSE = sqrt(sum(abs(est_ang - sig_angle).^2) / length(sig_angle));
end