function rmse = get_rmse(doa, est_doa)

    rmse = sqrt(sum(vec(min(abs(vec(doa) - vecH(est_doa)).').^2))/length(doa));
 
end