function sp = capon(recv, target_num, spatial_steer_mat)
    % ref: https://blog.csdn.net/pwang95/article/details/106343667/

    cov_mat = recv * recv'; 
    sp = 1./abs(sum((spatial_steer_mat' / cov_mat).' .* spatial_steer_mat));
end