function est_doa = root_music(recv, target_num,d)
    % ref: https://www.comm.utoronto.ca/~rsadve/Notes/DOA.pdf
    cov_mat = recv * recv';
    [eig_vector, eig_value] = eig(cov_mat);
    eig_value = abs(diag(eig_value));
    [~, sort_idx] = sort(eig_value, 'descend');
    noise_space = eig_vector(:, sort_idx(target_num + 1:end)); 
    C = noise_space * noise_space';
    coef = zeros(2 * size(C, 1) - 1, 1);

    for idx =- (size(C, 1) - 1):size(C, 1) - 1
        coef(idx + size(C, 1)) = sum(diag(C, idx));
    end

    a = roots(coef);
    a = a(abs(a) < 1);
    [~, sort_idx] = sort(abs(abs(a) - 1), 'ascend');
    z = a(sort_idx(1:target_num));
    est_doa = sort(-asind(angle(z) / (2 * pi * d)), 'ascend');
end
