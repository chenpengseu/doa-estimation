function sp = somp(recv, target_num, dic_mat)

    r_tmp = recv;
    idx_set = [];

    for idx = 1:target_num
        [~, max_idx] = max(sum(abs(r_tmp' * dic_mat).^2));
        idx_set = [idx_set; max_idx];
        r_tmp = recv - dic_mat(:, idx_set) * pinv(dic_mat(:, idx_set)) * recv;
    end

    sp_value = sum((abs(pinv(dic_mat(:, idx_set)) * recv).^2).');
    sp = zeros(size(dic_mat, 2), 1);
    sp(idx_set) = sp_value;
end
