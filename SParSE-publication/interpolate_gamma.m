function g_opt = interpolate_gamma(gammas_mat, counters, target, k, k_default, abs_tol, Kce, interp_th, dirName, retry)

fprintf('\n*** Starting interpolation function ***\n');
%% start interpolation
M           = length(gammas_mat(:,1));
markersz    = 8;
res_size    = 20;
tcounter    = target*Kce;
num_interp  = 4;

[counters, i_c] = sort(counters);
gammas_mat      = gammas_mat(:,i_c);
nz_ind          = find(counters>interp_th,1,'first');
if isempty(nz_ind)
    nz_ind      = find(counters> 0,1,'first');
end
nu_ind          = find(counters<(Kce-interp_th),1,'last');
if isempty(nu_ind)
    nu_ind      = find(counters< Kce, 1,'last');
end

if counters(nz_ind) > tcounter
    nz_ind = nz_ind-1;
end
if counters(nu_ind) < tcounter
    nu_ind = nu_ind+1;
end

counters   = counters(nz_ind:nu_ind);
gammas_mat = gammas_mat(:,nz_ind:nu_ind);

%% select viable interpolation points
dist       = abs(counters - tcounter);
[~, ind_d] = sort(dist);
if length(counters) > num_interp
    % divide between UP and OP
    u_ind  = find(counters<tcounter,1,'last');
    o_ind  = find(counters>tcounter,1,'first');
    
    % at least one on either side
    inds_temp = ind_d(1:num_interp);
    if min(inds_temp) > u_ind
        inds_temp(end) = u_ind;
    end
    if max(inds_temp) < o_ind
        inds_temp(end) = o_ind;
    end
    counters_d   = counters(inds_temp);
    gammas_mat_d = gammas_mat(:,inds_temp);
else
    % use all points
    counters_d   = counters(ind_d);
    gammas_mat_d = gammas_mat(:,ind_d);
end

%% plot the current gamma values and their counters
figure(1)
hold on
% plot computed gamma values from stages 1 and 2
plot(gammas_mat_d', counters_d','o','markersize',markersz)
% plot horizontal line at the target event
txs = linspace(min(min(gammas_mat_d)), max(max(gammas_mat_d)), res_size);
tys = linspace(tcounter, tcounter,res_size);
plot(txs,tys,':','linewidth',2.5,'Color',[0.9451 0.77255 0.07451]);
grid on

%% exponential fit with polyfit
counters_log   = log(counters_d);
poly_coiffs    = zeros(M, 2);
poly_coiffs_up = zeros(M, 2);
poly_coiffs_op = zeros(M, 2);
exps           = zeros(M, res_size);
xs             = zeros(M, res_size);
g_opt          = zeros(1,M);
% find farthest up and op
[counters_s, ind_s] = sort(counters_d);
gammas_mat_s        = gammas_mat_d(:,ind_s);
% UP
up_far = -1;
if sum(counters_s < tcounter) > 1 && length(counters_s)>2
    up_far = find(counters_s < tcounter);
    up_far = up_far(1);
    
    if up_far > 0
        gamma_mat_up_far = gammas_mat_s;
        gamma_mat_up_far(:,up_far) = [];
        counters_up_far  = counters_s;
        counters_up_far(up_far) = [];
        counters_log_upf = log(counters_up_far);
    end
end
% OP
op_far = -1;
if sum(counters_s > tcounter) > 1 && length(counters_s)>2
    op_far = find(counters_s > tcounter);
    op_far = op_far(end);
    
    if op_far>0
        gamma_mat_op_far = gammas_mat_s;
        gamma_mat_op_far(:,op_far) = [];
        counters_op_far  = counters_s;
        counters_op_far(op_far) = [];
        counters_log_opf = log(counters_op_far);
    end
end

% Fit polynomail and translate to exponential func
for j=1:M
    poly_coiffs(j,:) = polyfit(gammas_mat_d(j,:), counters_log, 1);
    
    % check for the goodness of fit by removing the farthest point
    % http://stats.stackexchange.com/questions/142248/difference-between-r-square-and-rmse-in-linear-regression\
    % R2 is the fraction of the total sum of squares that is 'explained by' the regression.
    y_hat = exp(poly_coiffs(j,2)) * exp(poly_coiffs(j,1)*gammas_mat_d(j,:));
    sse   = sum((counters_d - y_hat).^2);
    y_bar = sum(counters_d)/length(counters_d);
    tss   = sum((counters_d - y_bar).^2);
    r2    = 1 - sse/tss;
    fprintf('Rxn %d R2 goodness of fit: %f \n', j, r2);
    
    % remove farthest up and repeat
    if up_far > 0
        poly_coiffs_up(j,:) = polyfit(gamma_mat_up_far(j,:), counters_log_upf, 1);
        y_hat_up  = exp(poly_coiffs_up(j,2)) * exp(poly_coiffs_up(j,1)*gamma_mat_up_far(j,:));
        sse_up    = sum((counters_up_far - y_hat_up).^2);
        y_bar_up  = sum(counters_up_far)/length(counters_up_far);
        tss_up    = sum((counters_up_far - y_bar_up).^2);
        r2_up     = 1 - sse_up/tss_up;
        switch_up = r2_up > r2;
        fprintf('Rxn %d R2 goodness of fit with farthest up removal: %f \n', j, r2_up);
    end
    
    % remove farthest op and repeat
    if op_far > 0
        poly_coiffs_op(j,:) = polyfit(gamma_mat_op_far(j,:), counters_log_opf, 1);
        y_hat_op  = exp(poly_coiffs_op(j,2)) * exp(poly_coiffs_op(j,1)*gamma_mat_op_far(j,:));
        sse_op    = sum((counters_op_far - y_hat_op).^2);
        y_bar_op  = sum(counters_op_far)/length(counters_op_far);
        tss_op    = sum((counters_op_far - y_bar_op).^2);
        r2_op     = 1 - sse_op/tss_op;
        switch_op = r2_op > r2;
        fprintf('Rxn %d R2 goodness of fit with farthest op removal: %f \n', j, r2_op);
    end
    
    % which interpolation to choose
    % r2    = 1
    % r2_up = 2
    % r2_op = 3
    r_flag  = 1;
    switch sign(op_far)+sign(up_far)
        case -2
            % only 1 data for up and 1 data for op
        case 0
            % find op_far or up_far
            if up_far > op_far
                if switch_up
                    r_flag = 2;
                end
            else
                if switch_op
                    r_flag = 3;
                end
            end
        case 2
            % both op_far and up_far indices exist
            if switch_up
                if switch_op
                    if r2_up > r2_op
                        r_flag = 2;
                    else
                        r_flag = 3;
                    end
                else
                    r_flag = 2;
                end
            elseif switch_op
                r_flag = 3;
            end
        otherwise
            error('unknown case for deciding interpolant')
    end
    
    switch r_flag
        case 1
            g_opt(j)  = (log(tcounter)-poly_coiffs(j,2))/poly_coiffs(j,1);
            xs(j,:)   = linspace(min(gammas_mat_d(j,:)), max(gammas_mat_d(j,:)), res_size);
            exps(j,:) = exp(poly_coiffs(j,2)) * exp(poly_coiffs(j,1)*xs(j,:));
        case 2
            fprintf('Fit is better after removing farthest up. Computing new fit\n');
            g_opt(j)  = (log(tcounter)-poly_coiffs_up(j,2))/poly_coiffs_up(j,1);
            xs(j,:)   = linspace(min(gamma_mat_up_far(j,:)), max(gamma_mat_up_far(j,:)), res_size);
            exps(j,:) = exp(poly_coiffs_up(j,2)) * exp(poly_coiffs_up(j,1)*xs(j,:));
            plot(gammas_mat_d(:,ind_s(up_far)), counters_d(ind_s(up_far)),'kh','MarkerSize',markersz,...
                'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
        case 3
            fprintf('Fit is better after removing farthest op. Computing new fit\n');
            g_opt(j)  = (log(tcounter)-poly_coiffs_op(j,2))/poly_coiffs_op(j,1);
            xs(j,:)   = linspace(min(gamma_mat_op_far(j,:)), max(gamma_mat_op_far(j,:)), res_size);
            exps(j,:) = exp(poly_coiffs_op(j,2)) * exp(poly_coiffs_op(j,1)*xs(j,:));
            plot(gammas_mat_d(:,ind_s(op_far)), counters_d(ind_s(op_far)),'kh','MarkerSize',markersz,...
                'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
        otherwise
            error('unknown r_flag type for projecting g_opt')
    end
end

%% take care of interpolated values exceeding over/under perturbation
u_ind     = find(counters_s < tcounter, 1, 'last');
gammas_up = gammas_mat_s(:,u_ind);
o_ind     = find(counters_s > tcounter, 1, 'first');
gammas_op = gammas_mat_s(:,o_ind);
% reg_check = 1;
for j=1:M
    % compute weighted average with go_j
    dist    = abs([counters_s(o_ind) counters_s(u_ind)]-tcounter);
    dist    = dist / sum(dist); % normalize
    upflag  = 0;
    opflag  = 0;
    
    if g_opt(j) > 1
        if gammas_up(j) < gammas_op(j) % normal biasing
            if gammas_up(j) > g_opt(j)
                upflag = 1;
            elseif gammas_op(j) < g_opt(j)
                opflag = 1;
            end
        else % inverse biasing
            if gammas_up(j) < g_opt(j)
                upflag = 1;
            elseif gammas_op(j) > g_opt(j)
                opflag = 1;
            end
        end
    else
        if gammas_up(j) > gammas_op(j) % normal biasing
            if gammas_up(j) < g_opt(j)
                upflag = 1;
            elseif gammas_op(j) > g_opt(j)
                opflag = 1;
            end
        else % inverse biasing
            if gammas_up(j) > g_opt(j)
                upflag = 1;
            elseif gammas_op(j) < g_opt(j)
                opflag = 1;
            end
        end
    end
    %     end
    
    if upflag || opflag
        if upflag
            fprintf('UP ind %d: g_opt = %f \t g_up = %f \n', j, g_opt(j), gammas_up(j));
        else
            fprintf('OP ind %d: g_opt = %f \t g_op = %f \n', j, g_opt(j), gammas_op(j));
        end
        g_opt(j) = gammas_up(j) * dist(1) + gammas_op(j) * dist(2);
        fprintf('ind %d: new g_opt = %f \n', j, g_opt(j));
    end
end

%% plot the exponential fit and projection of the target event
plot(xs', exps', '--', 'LineWidth',1.5)
plot(g_opt, ones(1,M) * tcounter, '^', 'markersize', markersz,'Color',[0.847 0.1608 0])
hold off

% print ans save figure
k_factors = k_default./k;
for j=1:M
    fprintf(['projected optimal biasing parameter candidates for R',num2str(j),', wrt to k_default: ', num2str(g_opt(j),'%.5f '),'\n'])
    g_opt(j) = g_opt(j) * k_factors(j);
    fprintf(['candidate rates for the next iteration: ', num2str(g_opt(j)*k(j),'%.5f  '),'\n'])
end

saveas(gcf,[dirName,'\interpolation-target',num2str(target,'%.2f'),'_abstol',num2str(abs_tol),'_k',num2str(k_default,'_%.5f'),'_',num2str(retry'),'.jpg'])
saveas(gcf,[dirName,'\interpolation-target',num2str(target,'%.2f'),'_abstol',num2str(abs_tol),'_k',num2str(k_default,'_%.5f'),'_',num2str(retry'),'.fig'])
close gcf