function [gammas, flag_leap_type] = weighted_average(counters, gammas_mat, target_counter, Kce, k, k_default, interp_th)

gammas         = -1;
flag_leap_type = -1;
if length(counters)==1
    return
end


%% see if both UP and OP exist
counters2    = counters(counters>interp_th & counters<(Kce-interp_th));
gammas_mat2  = gammas_mat(:,counters>interp_th & counters<(Kce-interp_th));
up_valid_sum = sum(counters2 < target_counter);
op_valid_sum = sum(counters2 > target_counter);

%% both UP and OP are recorded
if up_valid_sum > 0 && op_valid_sum > 0
    [counters2_s, i_c2] = sort(counters2);
    gammas_mat2_s       = gammas_mat2(:,i_c2);

    % look at the distance to determine whether to run interpolation or one more leaping
    u_ind_interp = find(counters2_s < target_counter, 1, 'last');
    o_ind_interp = find(counters2_s > target_counter, 1, 'first');
    
    counters2_up = counters2_s(u_ind_interp);
    counters2_op = counters2_s(o_ind_interp);
    c2u_dist     = target_counter - counters2_up;
    c2o_dist     = counters2_op - target_counter;
    
    
    % proceed with last leap
    flag_leap_type = 1;
    u_weight       = c2o_dist/(c2o_dist+c2u_dist);
    o_weight       = c2u_dist/(c2o_dist+c2u_dist);
    % find weighted average
    gammas         = gammas_mat2_s(:,u_ind_interp)*u_weight + gammas_mat2_s(:,o_ind_interp)*o_weight;
    gammas         = gammas'.*k_default./k;
    return
end

%% both UP and OP are invalid -> bisection
if up_valid_sum == 0 && op_valid_sum == 0 && min(counters)<=interp_th && max(counters)>=(Kce-interp_th)
    % proceed leaping with weighted average
    flag_leap_type    = 2; 
    [counters_s, i_c] = sort(counters);
    u_ind_bisect      = find(counters_s < target_counter, 1, 'last');
    o_ind_bisect      = find(counters_s > target_counter, 1, 'first');
    % half of the two
    gammas            = (gammas_mat(:,i_c(u_ind_bisect)) + gammas_mat2_s(:,i_c(o_ind_bisect)))/2;
    gammas            = gammas'.*k_default./k;
    return
end

%% 2 or more points on one side and none on the other
% no valid OP
if up_valid_sum > 0 && op_valid_sum ==0 && max(counters)>=(Kce-interp_th)&& length(counters)>=3
    fprintf('\tObserved: two more more valid UP and all invalid OP\n')
    % look at the distance to pick points for weighted average
    [counters2_s, i_c2] = sort(counters2);
    gammas_mat2_s       = gammas_mat2(:,i_c2);
    u_ind_wave          = find(counters2_s < target_counter, 1, 'last');
    
    [counters_s, i_c]   = sort(counters);
    gammas_mat_s        = gammas_mat(:,i_c);
    o_ind_wave          = find(counters_s > target_counter, 1, 'first');
    
    counters_up         = counters2_s(u_ind_wave);
    counters_op         = counters_s(o_ind_wave);
    cu_dist             = target_counter - counters_up;
    co_dist             = counters_op - target_counter;
     
    % proceed leaping with weighted average
    flag_leap_type      = 3;
    u_weight            = co_dist/(co_dist+cu_dist);
    o_weight            = cu_dist/(co_dist+cu_dist);
    % find weighted average
    gammas              = gammas_mat2_s(:,u_ind_wave)*u_weight + gammas_mat_s(:,o_ind_wave)*o_weight;
    gammas              = gammas'.*k_default./k;
    return
end

% no valid UP
if op_valid_sum > 0 && up_valid_sum ==0 && min(counters)<=interp_th && length(counters)>=3
    fprintf('\tObserved: all invalid UP and two or more valid OP\n')
    % look at the distance to pick points for weighted average
    [counters_s, i_c]   = sort(counters);
    gammas_mat_s        = gammas_mat(:,i_c);
    u_ind_wave          = find(counters_s < target_counter, 1, 'last');
    
    [counters2_s, i_c2]  = sort(counters2);
    gammas_mat2_s       = gammas_mat2(:,i_c2);
    o_ind_wave          = find(counters2_s > target_counter, 1, 'first');
    
    counters_up         = counters_s(u_ind_wave);
    counters_op         = counters2_s(o_ind_wave);
    cu_dist             = target_counter - counters_up;
    co_dist             = counters_op - target_counter;
     
    % proceed leaping with weighted average
    flag_leap_type      = 3;
    u_weight            = co_dist/(co_dist+cu_dist);
    o_weight            = cu_dist/(co_dist+cu_dist);
    % find weighted average
    gammas              = gammas_mat_s(:,u_ind_wave)*u_weight + gammas_mat2_s(:,o_ind_wave)*o_weight;
    gammas              = gammas'.*k_default./k;    
end