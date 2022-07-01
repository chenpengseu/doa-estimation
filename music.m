function sp = music(recv, target_num, spatial_steer_mat)
    % ref: https://en.wikipedia.org/wiki/MUSIC_(algorithm)

    cov_mat = recv * recv';
    [eig_vector, eig_value] = eig(cov_mat);
    eig_value = abs(diag(eig_value)); 
    [~, sort_idx] = sort(eig_value, 'descend');
    noise_space = eig_vector(:, sort_idx(target_num+1:end)); 
    sp = 1./sum(abs(noise_space'*spatial_steer_mat).^2);
end
