function est_doa = ESPRIT(recv, target_num,d)
    % ref: https://en.wikipedia.org/wiki/Estimation_of_signal_parameters_via_rotational_invariance_techniques
    cov_mat = recv * recv';
    [U, E, V] = svd(cov_mat);
    sub_mat1 = U(1:size(cov_mat, 1) - 1, 1:target_num);
    sub_mat2 = U(2:size(cov_mat, 1), 1:target_num);
    P = sub_mat1 \ sub_mat2;
    est_doa = asind(angle(eig(P)) / (2 * pi * d));
end
