function doa = generate_doa(param)
    while (1) 
        if param.doa_min_spacing * param.K >= (param.doa_max - param.doa_min) * 0.7
            error('The targets are too close!');
        end 
 
        doa = sort(rand(param.K, 1) * (param.doa_max - param.doa_min) + param.doa_min, 'ascend');
        if param.K > 1
            if (min(abs(doa(2:end) - doa(1: end - 1))) >= param.doa_min_spacing) 
                break;
            end 
        else
            break;
        end 
    end
end