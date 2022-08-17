function recv = generate_recv(param, sig_len, doa, SNR)
    % signal and noise power
    Ps = ones(param.K, 1); % signal power 
    Pn = sum(Ps) / param.K / db2pow(SNR); 

    % generate the signals
    sig = bsxfun(@times, sqrt(Ps.' / 2), randn(sig_len, param.K) + 1j * randn(sig_len, param.K)); % target signals
    noise = sqrt(Pn / 2) * (randn(param.N, sig_len) + 1j * randn(param.N, sig_len));
    recv = get_steervec(param.N, param.d, deg2rad(doa)) * sig.' + noise;
end