clc
clear

% Min Roh
% Last updated 08/22/2016

%********** SIRS system description **********
% S + I -> 2I
% I     -> R
% R     -> S
%
% X0 = [S(t0) I(t0) R(t0)] = [100 1 0]
%*********************************************

% beta_min  = 0.005;
% beta_max  = 0.15;
% gamma_min = 0.5;
% gamma_max = 4.0;
% omega_min = 0.1;
% omega_max = 3;

% script for plotting 3D SSA output for SIRS
clear
clc
close all

%% final reaction rates
final_rates = zeros(30,4);

%% load mat file from 20x20 grid SSA simulations
load p_mat

%% choose a folder with simulation data
d_name = uigetdir('..','Choose a directory with SParSE++ simulation data');
curr_dir = pwd;
cd(d_name);
files = dir('*all.mat');
if isempty(files)
    error('Missing a file with ensemble simulation information')
elseif length(files)>1
    [filename, pathname] = uigetfile('*all.mat','Multiple files exist. Choose one desired output file.');
    load(fullfile(pathname, filename));
else
    load(files.name);
    [pathname,filename,ext] = fileparts(files.name);
    filename = strcat(filename, ext);
end
fprintf('Loading file %s\n',filename);
cd(curr_dir);

%% process sim information from the filename
sim_expression = '(\d.\d)';
m = regexp(filename,sim_expression,'match');
target = str2double(m{1});
abstol = str2double(m{2});
assert(isnumeric(target))
assert(isnumeric(abstol))

%% custom color arrangement
orange  = [1 0.5 0.2];
sc  = orange;
sce = orange;
ic  = 'w';
ice = orange;
icl = orange;
ec  = 'r';
ece = 'r';

%% figure config
alpha_factor = 0.05;
width        = 950;
height       = 600;
grid_size    = 50;
nf_rot       = 400;
klen         = 31;
kslen        = 30;

%% rotation and label config
az           = -133.5;
el           = 20;
view_angle1  = [az el];
xl_pos       = [-0.72 -3.95  1.55];
yl_pos       = [-0.73 -0.95  1.4];
zl_pos       = [-0.02  5.27  0.44];
tb_pos       = [0.56 0.8 0.22 0.1];


%% interpolate to a larger grid
load SSA_params.mat
betas_interp   = linspace(beta_min,  beta_max,  grid_size);
gammas_interp  = linspace(gamma_min, gamma_max, grid_size);
omegas_inpterp = linspace(omega_min, omega_max, grid_size);
[xi,yi,zi]     = meshgrid(betas_interp, gammas_interp, omegas_inpterp);
p_mat_interp   = interp3(betas,gammas,omegas,p_mat,xi,yi,zi);
beta_min_sim  = min(min(kp_mat(:,1)), beta_min);
beta_max_sim  = max(max(kp_mat(:,1)), beta_max);
gamma_min_sim = min(min(kp_mat(:,2)), gamma_min);
gamma_max_sim = max(max(kp_mat(:,2)), gamma_max);
omega_min_sim = min(min(kp_mat(:,3)), omega_min);
omega_max_sim = max(max(kp_mat(:,3)), omega_max);


%% final reaction counter
fr_counter = 0;

%% figure with finer mesh grid and 3d interpolation
hFig = figure(1);
set(hFig, 'Position', [300 250 width height])
axes1 = axes('Parent',hFig);
grid on
hold on

%% draw reference walls
[X, Y] = meshgrid(betas_interp, gammas_interp);
Z = zeros(grid_size, grid_size);
for i=1:grid_size
    Z(1:end) = omegas_inpterp(i);
    surf(Y,X,Z,p_mat_interp(:,:,i)','LineStyle','none')
    alpha(alpha_factor);
end

set(gca,'xdir','reverse')
view(axes1,view_angle1);
[xb, yg, zo] = meshgrid(gammas_interp,betas_interp, omegas_inpterp);
[f, v, c] = isosurface(xb, yg, zo, p_mat_interp, target, p_mat_interp);
patch('Vertices', v, 'Faces', f, 'FaceVertexCData', c, 'FaceColor','interp', 'edgecolor', 'interp');
alpha(alpha_factor);


%% customization
set(gca, 'clim', [0 1])
colorbar('location','EastOutside','Position',[.95 .112 .015 .815])
cb = findobj(gcf,'Type','axes','Tag','Colorbar');
cbIm = findobj(cb,'Type','image');
alpha(cbIm,alpha_factor+0.5)
axis tight

% shading interp
xlabel('\beta','FontSize',16)
ylabel('\gamma','FontSize',16)
zlabel('\omega','FontSize',16)

hxl = get(get(gca,'xlabel'),'Position');
hyl = get(get(gca,'ylabel'),'Position');
hzl = get(get(gca,'zlabel'),'Position');
set(get(gca,'xlabel'),'Position',xl_pos);
set(get(gca,'ylabel'),'Position',yl_pos);
set(get(gca,'zlabel'),'Position',zl_pos);

%% SParSePE output plotting
% create a cell array to hold individual trajectories
kpe_len      = length(kp_end);
kp_cell      = cell(1, kpe_len);
kp_cell_len  = zeros(1,kpe_len);
% exception in length calculation for the first element
kp_cell{1}     = kp_mat(1:kp_end(1),:);
kp_cell_len(1) = kp_end(1);
final_rates(1,:) = kp_mat(kp_end(1),:);
% rest of the cell entries
for i=2:kpe_len
    kp_cell_len(i) = kp_end(i) - kp_end(i-1);
    kp_cell{i}  = kp_mat((kp_end(i-1)+1):kp_end(i),:);
    final_rates(i,:) = kp_mat(kp_end(i),:);
end


% see if failed indices exist
failed_expression = '\w*.\d_abs\w.\d*';
m = regexp(filename,failed_expression,'match');
failed_filename = [char(m),'_failed.mat'];
fprintf('Checking for failure in convergence among chosen simulation ensemble\n')
fprintf('\t %s\n',failed_filename);
if exist(failed_filename,'file')
    fprintf('Failed simulation(s) found\n');
    load(failed_filename);
    if ~isempty(failed_ind) && failed_ind(1) == 0
        failed_ind = failed_ind(2:end);
    end
    if ~isempty(failed_ind)
        failed_ind
        fi_len = length(failed_ind);
        failed_cell_ind = zeros(fi_len,1);
        for fi = 1:fi_len
            failed_cell_ind(fi) = find(kp_end>failed_ind(fi),1,'first');
        end
        
        kp_cell(failed_cell_ind)     = [];
        kp_cell_len(failed_cell_ind) = [];
    end
else
    fprintf('No failed simulation(s) found\n');
end
% readjust parameters
kp_inds     = 1:length(kp_cell);
kpl_max     = max(kp_cell_len);

%% plot first points
% find ones with length == 1
kp_1s_inds = find(kp_cell_len == 1);
if ~isempty(kp_1s_inds)
    kp1i_len   = length(kp_1s_inds);
    kpmat_1s   = zeros(kp1i_len, 4);
    for i=1:kp1i_len
        kpmat_1s(i,:) = kp_cell{kp_1s_inds(i)};
    end
    kp_1s_beta   = kpmat_1s(:,2);
    kp_1s_gamma  = kpmat_1s(:,1);
    kp_1s_omega  = kpmat_1s(:,3);
    kp_1s_p      = kpmat_1s(:,4);
    
    % delete identified data
    kp_cell(kp_1s_inds)     = [];
    kp_inds = 1:length(kp_cell);
    kp_cell_len(kp_1s_inds) = [];
    
    % separate into target and not for kp_1s
    t_1s = find(abs(kp_1s_p-target)<abstol);
    kp1_inds = setdiff(1:kp1i_len, t_1s);
    plot3(kp_1s_beta(kp1_inds),kp_1s_gamma(kp1_inds),kp_1s_omega(kp1_inds),'s','MarkerFaceColor',sc,'MarkerEdgeColor',sce)
    if ~isempty(t_1s)
        plot3(kp_1s_beta(t_1s),kp_1s_gamma(t_1s), kp_1s_omega(t_1s),'s','MarkerFaceColor',ec,'MarkerEdgeColor',ece)
    end
end
kp_inds_len = length(kp_inds);

% plot title and rest of first points
title_str = ['Target: ',num2str(target*100),'%  AbsTol: ', num2str(abstol*100),'%  Init # Param: ',num2str(kslen),'  Tot # Param: ', num2str(kp_end(end))];
th= title(title_str,'FontSize',14);
for i=1:kp_inds_len
    kp_cell_i = kp_cell{i};
    plot3(kp_cell_i(1,2),kp_cell_i(1,1),kp_cell_i(1,3),'s','MarkerFaceColor',sc,'MarkerEdgeColor',sce)
end
pause(.05)

% intermediate points
counter = 2;
while counter <= kpl_max
    kp_ploti = zeros(1,4);
    kp_plot_prev = zeros(1,4);
    kp_ploti_end = zeros(1,4);
    kp_plot_prev_end = zeros(1,4);
    
    for i=1:kp_inds_len
        if kp_cell_len(i)>counter
            kp_cell_i = kp_cell{i};
            kp_ploti(end+1,:) = kp_cell_i(counter,[2 1 3 4]);
            kp_plot_prev(end+1,:) = kp_cell_i(counter-1,[2 1 3 4]);
        elseif kp_cell_len(i)==counter
            kp_cell_i = kp_cell{i};
            kp_ploti_end(end+1,:) = kp_cell_i(counter,[2 1 3 4]);
            kp_plot_prev_end(end+1,:) = kp_cell_i(counter-1,[2 1 3 4]);
        end
    end
    % get rid of the first row of 0s
    kp_ploti         = kp_ploti(2:end,:);
    kp_plot_prev     = kp_plot_prev(2:end,:);
    kp_ploti_end     = kp_ploti_end(2:end,:);
    kp_plot_prev_end = kp_plot_prev_end(2:end,:);
    % intermediate points
    for j=1:length(kp_ploti(:,1))
        plot3([kp_plot_prev(j,1);kp_ploti(j,1)],[kp_plot_prev(j,2);kp_ploti(j,2)],[kp_plot_prev(j,3);kp_ploti(j,3)],'--','Color',icl);
        plot3(kp_ploti(j,1),kp_ploti(j,2),kp_ploti(j,3),'s','MarkerFaceColor',ic,'MarkerEdgeColor',ice);
    end
    
    % end points that reached the target
    for j=1:length(kp_ploti_end(:,1))
        plot3([kp_plot_prev_end(j,1);kp_ploti_end(j,1)],[kp_plot_prev_end(j,2);kp_ploti_end(j,2)],[kp_plot_prev_end(j,3);kp_ploti_end(j,3)],'--','Color',icl);
        plot3(kp_ploti_end(j,1),kp_ploti_end(j,2),kp_ploti_end(j,3),'s','MarkerFaceColor',ec,'MarkerEdgeColor',ece);
    end
    counter = counter+1;
    pause(0.05)
end
% last points
inds_last    = find(kp_cell_len==kpl_max);
kp_ploti     = zeros(1,4);
kp_plot_prev = zeros(1,4);
for i=1:length(inds_last)
    kp_cell_i = kp_cell{inds_last(i)};
    kp_ploti(end+1,:) = kp_cell_i(kpl_max,[2 1 3 4]);
    kp_plot_prev(end+1,:) = kp_cell_i(kpl_max-1,[2 1 3 4]);
end
kp_ploti     = kp_ploti(2:end,1:4);
kp_plot_prev = kp_plot_prev(2:end,1:4);
for j=1:length(kp_ploti(:,1))
    plot3([kp_plot_prev(j,1);kp_ploti(j,1)],[kp_plot_prev(j,2);kp_ploti(j,2)],[kp_plot_prev(j,3);kp_ploti(j,3)],'--','Color',orange);
    plot3(kp_ploti(j,1),kp_ploti(j,2),kp_ploti(j,3),'rs','MarkerFaceColor','r');
end

set(th, 'Visible','off');
nf = 400;
for i=1:nf
    view(axes1, [az+(i/nf)*180 el]);
    pause(0.05)
end
hold off
