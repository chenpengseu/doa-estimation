function steer_vec = get_steervec(ant_num, ant_dist, spatial_angle)
%get_steervec - Get the steerign vector 
% exp(1j*2*pi*nd/lambda*sin(theta))

steer_vec = exp(1j*2*pi*[0:ant_num-1].'*ant_dist*sin(spatial_angle(:).'));
    
end