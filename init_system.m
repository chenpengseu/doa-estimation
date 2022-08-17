function param = init_system()
    param.N = 16; % antenna number
    param.d = 0.5; % normalized antenna distance
    param.K = 3; % target number
    param.doa_min_spacing = 100 / param.N*2; % deg, the minimal target spacing
    param.doa_min = -60; % deg, the detection angle
    param.doa_max = 60; % deg, the detection angle
end