% test speeding up SIRS
% start with the worst case
clc, clear, close all

% set random number stream
seed = 99901;
s    = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

% load the model file and overwrite target if desired
ofile_prefix = 'birthDeath';
switch ofile_prefix
    case 'revIsom'
        ic.x0    = [100 0];
        ic.t0    = 0;
        ic.tf    = 10;
        target_e = 30;
        
        fx0   = ic.x0(2);
        M     = 2;
        
        k1min = 0.1;
        k1max = 0.3;
        k2min = 0.3;
        k2max = 1;
    case 'birthDeath'
        ic.x0    = 40;
        ic.t0    = 0;
        ic.tf    = 100;
        target_e = 80;
        
        fx0   = ic.x0;
        M     = 2;
        
        k1min = 1;
        k1max = 1.7;
        k2min = 0.0125;
        k2max = 0.025;
        
    case 'SIRS'
        ic.x0    = [100 1 0];
        ic.t0    = 0;
        ic.tf    = 30;
        target_e = 50;
        
        fx0       = ic.x0(2);
        M         = 3;
        
        beta_min  = 0.005;
        beta_max  = 0.15;
        gamma_min = 0.5;
        gamma_max = 4.0;
        omega_min = 0.1;
        omega_max = 3;
end


%% simulation params
targets       = [0.4, 0.6, 0.8]; 
abs_thrhs     = [0.01 0.05 0.1];
simCounterLim = 20;
Kce           = 5e4; 
kslen         = 30;

%% play with this number for robustness/accuracy
% max_iter + 1 (initial CE) iterations will be run
max_iter      = 10;
% minimum number of CE stages required to decide leaping
min_CEstgs    = 2;
% number of stages (first row) and steps (second row) for leaping
size_leaps    = [50 20 5;
    7  5  3];

%% determine the direction of bias
if (fx0 > target_e)
    event_type = -1;
else
    event_type = 1;
end

%% start simulation
fprintf(['SParSE v2 simulation with ', ofile_prefix,' system\n']);

for ts = 1:length(targets)
    % load target
    target    = targets(ts);
    if target < 0.5
        % set interpolation threshold
        interp_th = target*Kce*0.05;
        % set leaping threshold
        leap_th   = target*Kce*0.01;
    else        
        % set interpolation threshold
        interp_th = (1-target)*Kce*0.05;
        % set leaping threshold
        leap_th   = (1-target)*Kce*0.01;
    end
    
    %% generate random sets of parameters as starting point
    switch ofile_prefix
        case {'revIsom','birthDeath'}
            ks = generateKs(kslen, k1min, k1max, k2min, k2max);
        case 'SIRS'
            ks = generateKs(kslen, beta_min, beta_max, gamma_min, gamma_max, omega_min, omega_max);
    end
    
    % set up the location and name of output file and diary
    date_str  = datestr(now,'mm-dd-hh');
    directory = ['.\Output\',ofile_prefix,'\', date_str, '_target',num2str(target),'_Kce',num2str(Kce),'_kslen',num2str(kslen)];
    mkdir(directory);
    
    for as = 1:length(abs_thrhs)
        % absolute error threshold
        abs_thrh = abs_thrhs(as);
        fprintf('Target: %g \tInterpolation Threshold: %g \tLeaping Threshold: %g \tAbs tol: %g\n', target, interp_th, leap_th,abs_thrh);
        
        % for storing results
        output_filename_base = [directory,'\',ofile_prefix,'_target',num2str(target),'_abs',num2str(abs_thrh)];
        
        failed_ind = -1;
        kp_mat     = zeros(1,M+1);
        kp_end     = 0;
        % record number of stages
        num_stages = zeros(kslen, 5); %first second op leap interpolation
        
        %% leaping parameters
        % record when leaping occured
        leap_indices  = [];
        
        for kc = 1:kslen
            %% start measuring total time
            tt_start = tic;
            
            %% set random number stream
            seed = seed + 100;
            s    = RandStream('mt19937ar','Seed',seed);
            RandStream.setGlobalStream(s);
   
            %% simulation parameter values
            simCounter = 1;
            success    = 0;
            op_flag    = 0;
            sim_flag   = 1;
            int_flag   = 0;
            k_default  = ks(kc,:);
            k          = k_default;
            k_old      = k;
            d_thrh     = min(target*0.8, 0.2);
            gammas     = ones(1,M);
            t_counter  = target * Kce;
            
            %% leaping parameters
            % flag to indicate whether to run normal CE or leap
            flag_CE_leap   = 0;
            % number of normal CE stages run since the last leap
            num_CE         = 0;
            % local start index for leaping
            leap_start_ind = 1;
            % total number of leaping
            leap_no        = 0;
            % used to indicate last leaping prior to interpolation stage
            flag_last_leap = 0;
            
            %% gamma interpolation parameters
            counters   = [];
            gammas_mat = zeros(M, 1);
            
            
            %% set up the location and name of output file and diary
            output_filename = [output_filename_base,'_k',num2str(k,'_%.5f'),'_seed',num2str(seed)];
            diary off
            diary([output_filename,'.txt'])
            
            %% start the simulation
            fprintf(['\n##### Sim ',num2str(kc),'\tInitial state: ', num2str(ic.x0),' with rates: ',num2str(k),' #####']);
            while(simCounter < simCounterLim)
                % non-0 and non-1 p^
                flag_valid_leap = 0;
                g_len           = length(gammas(:,1));
                for g = 1:g_len
                    %% first stage
                    t1                  = tic;
                    gammai              = gammas(g,:);
                    if flag_last_leap
                        num_stages(kc,4)    = num_stages(kc,4)+1;
                        leap_indices(end+1) = length(counters)+1;
                        fprintf(['\nfirst stage, last leaping using gamma: ', num2str(gammai),'\n']);
                    elseif flag_CE_leap
                        num_stages(kc,4)    = num_stages(kc,4)+1;
                        leap_indices(end+1) = length(counters)+1;
                        fprintf(['\nfirst stage with leaping using gamma: ', num2str(gammai),'\n']);
                    else
                        num_stages(kc,1) = num_stages(kc,1)+1;
                        fprintf(['\nfirst mCE stage index: ',num2str(g),' using gamma: ', num2str(gammai),'\n']);
                    end
                    [counter, localExts] = firstStage(ic, k, gammai, target_e, event_type, Kce, ofile_prefix);
                    t1e                  = toc(t1);
                    fprintf('first stage of CE took %.2f sec\n', t1e);
                    
                    % process and store result
                    distance        = target - counter/Kce;
                    kp_mat(end+1,:) = [k.*gammai counter/Kce];
                    counters(end+1) = counter;
                    fprintf(['Target event was reached ', num2str(counter), ' times with gamma ', num2str(gammai), '\n']);
                    fprintf('Distance between the target and current fraction is %.5f and p_hat is: %.5f\n', distance, counter/Kce);
                    
                    %% check for target probability
                    if abs(distance) < abs_thrh
                        simCounter = simCounterLim+1;
                        success    = 1;
                        sim_flag   = 0;
                        tt_end     = toc(tt_start);
                        k_final    = k.*gammai;
                        
                        kp_end(end+1) = length(kp_mat(:,1));
                        gammas_mat(:,length(counters)) = k_final./k_default;
                        
                        fprintf('\n\n*****Success: less than %0.3f distance remaining from target*****\n', abs_thrh);
                        fprintf(['final rate is ', num2str(k_final),'\n']);
                        fprintf(['total simulation time is ', num2str(tt_end),'\n\n']);
                        
                        save([output_filename,'.mat']);
                        diary off
                        break;
                    end
                    
                    %% record interpolation entry
                    % old code had num2str(length(counters))+1. does not affect code result but counter was off by 1 before
                    gammas_mat(:,length(counters)) = gammai.*k./k_default;
                    fprintf(['counter at ', num2str(length(counters)),': ', num2str(counter), '\n']);
                    fprintf(['gammas at ', num2str(length(counters)),' with respect to k_default: ', num2str(gammai.*k./k_default, '%.3f '), '\n']);
                    
                    if flag_last_leap
                        flag_last_leap = 0;
                        int_flag       = 1;
                        simCounter     = simCounterLim+1;
                        break;
                    end
                    
                    %% validation for leaping
                    % tester for leaping, once enough CE stages are run
                    if flag_CE_leap
                        flag_CE_leap = 0;
                        [gammas, flag_leap_type] = weighted_average(counters, gammas_mat, tcounter, Kce, k, k_default, interp_th);
                        
                        switch(flag_leap_type)
                            case 1 % run one last time
                                flag_last_leap = 1;
                                break;
                            case 2 % bisection
                                flag_CE_leap   = 1;
                                fprintf('\tEntering leaping with bisection\n');
                                break;
                            case 3 % weighted average
                                flag_CE_leap   = 1;
                                fprintf('\tEntering leaping with weighted average\n');
                                break;
                        end
                    else
                        leap_condition = 0;
                        tcounter       = target * Kce;
                        % check for isocline crossing
                        [gammas_check, flag_leap_type] = weighted_average(counters, gammas_mat, tcounter, Kce, k, k_default, interp_th);
                        
                        if gammas_check > 0 %anything but 0
                            gammas = gammas_check;
                        end
                        switch(flag_leap_type)
                            case 1 % run one last time
                                flag_last_leap = 1;
                                break;
                            case 2 % bisection
                                flag_CE_leap   = 1;
                                fprintf('\tEntering leaping with bisection\n');
                                break;
                            case 3 % weighted average
                                flag_CE_leap   = 1;
                                fprintf('\tEntering leaping with weighted average\n');
                                break;
                        end
                        
                        if ~flag_valid_leap && counter > leap_th && counter < (Kce-leap_th)
                            flag_valid_leap = 1;
                            num_CE          = num_CE + 1;
                            leap_condition  = 0;
                            if num_CE == min_CEstgs
                                leap_file_path = [directory,'\CEleap-target',num2str(target,'%.2f'),'_abstol',num2str(abs_thrh),'_k',num2str(k_default,'_%.5f'),'_',num2str(leap_no')];
                                if leap_start_ind>1
                                    [leap_condition, gammas_check] = checkLeap(target, k, k_default, Kce, size_leaps, leap_th, counters(leap_start_ind-1:end), gammas_mat(:,leap_start_ind-1:end), leap_file_path);
                                else
                                    [leap_condition, gammas_check] = checkLeap(target, k, k_default, Kce, size_leaps, leap_th, counters(leap_start_ind:end), gammas_mat(:,leap_start_ind:end), leap_file_path);
                                end
                                if leap_condition == 1
                                    gammas       = gammas_check;
                                    flag_CE_leap = 1;
                                    leap_no      = leap_no+1;
                                    break;
                                end
                            end
                        end
                    end
                    
                    %% normal CE stage - move on or compute IEs
                    if ~op_flag                 
                        [counters_s, ~] = sort(counters);
                        counters2  = counters_s(counters_s>interp_th);
                        counters2  = counters2(counters2<(Kce-interp_th));
                        if sign(distance) ~= event_type
                            fprintf('\n***** Change in bias direction: OP from UP ******\n');
                            if g == g_len
                                if sum(counters2 < target*Kce) == 0  % no valid UP
                                    op_flag = 1;
                                elseif sum(counters2 < target*Kce) > 0 && sum(counters2 > target*Kce) > 0% both UP and OP recorded
                                    fprintf('***** Exhausted gammas with both UP and OP records *****\nProceed with gamma interpolation\n');
                                    int_flag = 1;
                                    simCounter = simCounterLim + 1;
                                    break;
                                end
                            else
                                fprintf(['**Trying next gamma: ', num2str(gammas(g+1,:)), '**\n']);
                                continue
                            end
                        else
                            if sum(counters2 < target*Kce) > 0 && sum(counters2 > target*Kce) > 0% both UP and OP recorded
                                int_flag = 1;
                                fprintf('***** Second isocline crossing *****\nProceed with gamma interpolation\n');
                                simCounter = simCounterLim + 1;
                                break
                            else
                                k_old     = k;
                                gamma_old = gammai;
                                k         = k.*gammai;
                                gammas    = ones(1,M);
                                fprintf(['**Updating rate from ', num2str(k_old), ' to ', num2str(k),' and resetting gammas to 1s... **\n']);
                                
                                [ies, prcs, rhois] = findIE(localExts, Kce, distance, counter, event_type);
                                for l=1:length(ies)
                                    fprintf('IEs for rhos %.3f are %d with percs %.2f%%\n', rhois(l), ies(l), prcs(l)/Kce*100);
                                end
                                break
                            end
                        end
                    end
                    
                    
                    %% op flag
                    if abs(op_flag)                     
                        [counters_s, ~] = sort(counters);
                        counters2       = counters_s(counters_s>interp_th);
                        counters2       = counters2(counters2<(Kce-interp_th));
                        
                        if sign(distance) == event_type
                            fprintf('\n***** Change in bias direction: UP from OP ******\n');
                            if g == g_len
                                if sum(counters2 < target*Kce) > 0 && sum(counters2 > target*Kce) > 0% both UP and OP recorded
                                    fprintf('***** Exhausted gammas with both UP and OP records *****\nProceed with gamma interpolation\n');
                                    int_flag = 1;
                                    simCounter = simCounterLim + 1;
                                    break;
                                elseif sum(counters2 < target*Kce) == 0 % no valid UP
                                    op_flag      = 0;
                                    fprintf('::UP from OP:: without valid UP. Start leaping. \n');
                                    [gammas, ~]  = weighted_average(counters, gammas_mat, target*Kce, Kce, k, k_default, interp_th);
                                    flag_CE_leap = 1;
                                    break;
                                end % else no valid OP
                            else
                                fprintf(['**Trying next gamma: ', num2str(gammas(g+1,:)), '**\n']);
                                continue
                            end
                        else
                            if sum(counters2 < target*Kce) > 0 && sum(counters2 > target*Kce) > 0% both UP and OP recorded
                                int_flag = 1;
                                fprintf('***** Second isocline crossing *****\nProceed with gamma interpolation\n');
                                simCounter = simCounterLim + 1;
                                break
                            else
                                k_old  = k;
                                k      = k.*gammai;
                                gammas = ones(1,M);
                                fprintf(['::OP:: Updating rate from ', num2str(k_old), ' to ', num2str(k),' and resetting gammas to 1s... **\n']);
                            end
                        end
                        
                        % OP stage
                        t1_op   = tic;
                        fprintf(['::OP:: stage using gamma: ', num2str(gammas),' and k: ',num2str(k),'\n']);
                        localExts = opStage(ic, k, event_type, Kce, ofile_prefix);
                        t1e_op = toc(t1_op);
                        fprintf('::OP:: stage of CE took %.2f sec\n', t1e_op);
                        num_stages(kc,3) = num_stages(kc,3)+1;
                        
                        % determine OP IEs
                        [ies, prcs, rhois] = findIE(localExts, Kce, distance, counter, event_type);
                        for l=1:length(ies)
                            fprintf('::OP:: IEs for rhos %.3f are %d with percs %.2f%%\n', rhois(l), ies(l), prcs(l)/Kce*100);
                        end
                        break
                    end
                end
                
                % decide whether to leap or continue with CE
                if flag_last_leap && sim_flag
                    fprintf('### Both UP and OP are reached with last stage as leaping. Last leaping initiated with gamma from weighted average ###\n');
                    leap_start_ind = length(counters)+1;
                elseif flag_CE_leap && sim_flag
                    fprintf('\t### Leaping initiated, skip the second stage ###\n');
                    % reset number of CE
                    num_CE         = 0;
                    % reset leaping
                    leap_start_ind = length(counters)+1;
                elseif leap_condition == -1
                    % reset number of CE
                    num_CE         = 0;
                    % reset leaping
                    leap_start_ind = length(counters)+1;
                end
                
                %% second sim
                if (simCounter < simCounterLim) && ~flag_CE_leap && ~flag_last_leap
                    t2 = tic;
                    fprintf(['\tsecond stage using gamma: ', num2str(gammas),' and k: ',num2str(k),'\n']);
                    % COUNTER NOT NEEDED HERE
                    [num, denom, ~] = secondStage(ic, k, ies, Kce, ofile_prefix);
                    t2e = toc(t2);
                    fprintf('\tsecond stage of CE took %.2f sec\n', t2e);
                    num_stages(kc,2) = num_stages(kc,2)+1;
                    % update data and print info
                    if sign(distance) == event_type % equiv to op_flag = 0
                        gammas = num./denom;
                        if sum(sum(isnan(gammas)))>0 || sum(sum(gammas==0))
                            display('\tresetting 0 and nan entries to double.min')
                            gammas(isnan(gammas)) = realmin;
                            gammas(gammas==0) = realmin;
                        end
                    else
                        op_flag        = -1;
                        gammas = denom./num;
                        if sum(sum(isnan(gammas)))>0 || sum(sum(gammas==0)) || sum(sum(gammas > 20))>0
                            display('\tresetting 0 and nan entries to double.min')
                            gammas(isnan(gammas)) = realmin;
                            gammas(gammas==0) = realmin;
                            
                            display('\tresetting gamma>20 to 20 if applicable')
                            gammas(gammas>20) = 20;
                        end
                    end
                    for i=1:length(ies)
                        fprintf(['\tGamma for IE = ', num2str(ies(i)), ' is ', num2str(gammas(i,:)), '\n']);
                    end
                end
                
                simCounter = simCounter + 1;
            end
            
            if simCounter == simCounterLim
                fprintf('\n*** Sim Counter Limit has been reached. Add sim index to the list of failed indices and move on... ***\n\n');
                if failed_ind(1) == -1
                    failed_ind(1) = kc;
                else
                    failed_ind(end+1)= kc;
                end
                kp_end(end+1) = length(kp_mat(:,1));
                continue;
            end
            
            %% interpolation stage
            retry     = 1;
            retry_lim = 10;
            if sim_flag
                while retry <= retry_lim;
                    % interpolate and project optimal gamma estimate
                    fprintf(['\nInterpolation index: ',num2str(retry),': Compute candidate gammas from interpolation\n']);
                    tg               = tic;
                    gammas_opt       = interpolate_gamma(gammas_mat, counters, target, k, k_default, abs_thrh, Kce, interp_th, directory, retry);
                    fprintf(['Starting ',num2str(retry),' simulation from interpolated gammas: ',num2str(gammas_opt,'%.3f  '),'\n']);
                    % first stage
                    [counter, x_ext] = firstStage(ic, k, gammas_opt, target_e, event_type, Kce, ofile_prefix);
                    tge              = toc(tg);
                    fprintf('Interpolation stage took %.2f sec\n', tge);
                    num_stages(kc,5) = num_stages(kc,5)+1;
                    
                    % check for convergence
                    distance = target - counter/Kce;
                    fprintf(['\t## Target event was reached ', num2str(counter), ' times with gamma ', num2str(gammas_opt), '##\n']);
                    fprintf('\t## Distance between the target and current fraction is %.3f##\n', distance);
                    
                    counters(end+1)                = counter;
                    gammas_mat(:,length(counters)) = gammas_opt .* k./k_default;
                    kp_mat(end+1,:)                = [k.*gammas_opt counter/Kce];
                    
                    if abs(distance) < abs_thrh
                        sim_flag      = 0;
                        tt_end        = toc(tt_start);
                        k_final       = k.*gammas_opt;
                        kp_end(end+1) = length(kp_mat(:,1));
                        
                        fprintf('\n\t*****Success: less than %0.3f distance remaining from target*****\n\t     Exisiting Interpolation...\n', abs_thrh);
                        fprintf(['\t     Final rate is ', num2str(k_final),'\n']);
                        fprintf(['\t     Total simulation time is ', num2str(tt_end),'\n\n']);
                        
                        save([output_filename,'.mat']);
                        diary off
                        break;
                    else
                        retry = retry + 1;
                    end
                end
            end
            if retry > retry_lim
                fprintf('\n*** Interpolation Counter Limit has been reached. Add sim index to the list of failed indices and move on... ***\n\n');
                if failed_ind(1) == -1
                    failed_ind(1) = kc;
                else
                    failed_ind(end+1)= kc;
                end
                kp_end(end+1) = length(kp_mat(:,1));
            end
        end
        
        kp_mat = kp_mat(2:end,:);
        kp_end = kp_end(2:end);
        kp_end = kp_end - 1;
        output_filename_mat = [directory,'\',ofile_prefix,'_target',num2str(target),'_abs',num2str(abs_thrh)];
        if failed_ind(1) ~= -1
            failed_ind = failed_ind(2:end);
            save([output_filename_mat,'_failed.mat'],'failed_ind');
        end
        save([output_filename_mat,'_kslen',num2str(kslen),'_all.mat']);
    end
end