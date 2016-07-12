function [leap_condition, gammas_opt] = checkLeap(target, k, k_default, Kce, size_leaps, leap_th, counters, gammas_mat, file_path)
leap_condition = -1; % not qualified for leaping
gammas_opt     = -1; % not qualified for leaping

inds = intersect(find(counters>leap_th), find(counters<(Kce-leap_th)));
if length(inds) < 2
    error('not enough (<2) valid data is present\n')
end

counters   = counters(inds);
gammas_mat = gammas_mat(:, inds);
steps      = abs(diff(counters));
mean_step  = mean(steps);
[M,~]      = size(gammas_mat);

% number of steps CE will take, linear approximation
% use counter closest to target
counter_diff         = target*Kce - counters;
[min_cdiff, mcd_ind] = min(abs(counter_diff));
future_step          = floor(min_cdiff/mean_step);

% decide how much to leap
leap_ind             = find(future_step<=size_leaps(1,:),1,'last');
if leap_ind == length(size_leaps(1,:)) % not enough steps remain
    fprintf('\t### NOT qualified for leaping: mean_step:%g\tapprx_future_steps:%d\tcutoff step size: %d\n', mean_step,future_step, size_leaps(1,end));
    return
elseif isempty(leap_ind) % maximum leap
    fprintf('\t### Qualified for leaping: mean_step:%g\tapprx_future_steps:%d\tmultiplier: %d\n', mean_step,future_step, size_leaps(2,1));
    leap_multiplier = size_leaps(2,1);
else
    fprintf('\t### Qualified for leaping: mean_step:%g\tapprx_future_steps:%d\tmultiplier: %d\n', mean_step,future_step, size_leaps(2,leap_ind+1));
    leap_multiplier = size_leaps(2,leap_ind+1);
end
% turn on leaping
leap_condition = 1; 
% compute target probability
target_c       = counters(mcd_ind) + sign(counter_diff(mcd_ind))*mean_step * leap_multiplier;

if sign(counter_diff(mcd_ind))==1
    % UP
    if target_c > target*Kce
        target_c = target*Kce;
    end
else
    % OP
    if target_c < target*Kce
        target_c = target*Kce;
    end
end

% based on target_p and ks_mat, extrapolate next gamma using exponential
% function: counters = exp(a*ks+b) => log(counters)=a*gamma+b
counters_log = log(counters);
poly_coiffs  = zeros(M,2);
% for plotting
exps           = zeros(M, 20);
xs             = zeros(M, 20);
for j=1:M
    poly_coiffs(j,:) = polyfit(gammas_mat(j,:), counters_log, 1);
    gammas_opt(j)    = (log(target_c)-poly_coiffs(j,2))/poly_coiffs(j,1);
    % for plotting
    xs(j,:)          = linspace(min(gammas_opt(j),min(gammas_mat(j,:))), max(gammas_opt(j),max(gammas_mat(j,:))), 20);
    exps(j,:)        = exp(poly_coiffs(j,2)) * exp(poly_coiffs(j,1)*xs(j,:));
end
%% plot the current gamma values and their counters
figure(1)
hold on
% plot computed gamma values from stages 1 and 2
h_marker = plot(gammas_mat', counters','o','markersize',8);
% plot horizontal line at the target event
txs = linspace(min(min(min(gammas_mat)),min(min(gammas_opt))), max(max(max(gammas_mat)),max(max(gammas_opt))), 20);
tys = linspace(target_c, target_c, 20);
plot(txs,tys,':','linewidth',2.5,'Color',[0.9451 0.77255 0.07451]);
grid on
plot(gammas_opt, ones(1,M) * target_c, '^', 'markersize', 8,'Color',[0.847 0.1608 0])
plot(xs', exps', '--','LineWidth',1.5)
hold off
%% rescale with respect to the current k
gammas_opt = gammas_opt .* k_default./k;

%% save figure file
saveas(gcf,[file_path,'.jpg'])
saveas(gcf,[file_path,'.fig'])
close gcf
end
