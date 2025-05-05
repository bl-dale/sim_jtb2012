% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% Figure 6
% Plot pHs for different UL widths
% This is Figure 6 of the paper

function [] = create_fig6_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,pH_outs,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('Number of data file: ');
%n = length(all_data)

%strcat(sim_dir,sim_filename_base,'.mat')
%load(strcat(sim_dir,'/',sim_filename_base,'.mat'));

shades = [0    0.5000         0;  % green
         0         0    1.0000;   % blue
    1.0000         0         0;   % red
         0    0.7500    0.7500;   % cyan
    0.7500         0    0.7500;   % magenta
    162/255  120/255        0;   % gold
    0.2500    0.2500    0.2500];  % black

set(0,'DefaultAxesColorOrder',shades);
set(0,'defaultaxesfontsize',18);

figure;


for i=1:n   
    %load(strcat(sim_dir,'\',all_data(i).name));
    time   = cell2mat(times(i));
    X      = cell2mat(Xs(i));
    n_in   = cell2mat(n_ins(i));
    N      = cell2mat(Ns(i));
    n_buff = cell2mat(n_buffs(i));
    R      = cell2mat(Rs(i));
    R_inf  = cell2mat(R_infs(i));
    pH_out = cell2mat(pH_outs(i));

    d_EUS = (R_inf - R)*1e4;
    
    %label
    name = 'd = ' + string(round(d_EUS)) + ' \mum';
    
    n1 = (1+n_buff)*N + n_in; % one shell below membrane; n1+1 = @membrane
    tf_i = 100;
    tf_s = 100;
    depth = 50;      % depth of electrode inside (in microns)
    depth = 1e-4*depth;  % d in centimeters
    rad_in = (R/n_in)*[0:n_in];
    ind_electrode = find(rad_in >= R-depth,1); % inside

    time_aux = [-100];
    time_plot = [time_aux'; time];

    pHi = 3-log10(X(:,n1-(n_in-ind_electrode)));
    pHi_aux = [pHi(1)];
    pHi_plot = [pHi_aux'; pHi];

    pHs = 3-log10(X(:,n1+2));
    pHs_aux = [pHs(1)];
    pHs_plot = [pHs_aux'; pHs];
   
    for j = 1:size(time,1)-1
        time_minus = time(j);
        time_plus = time(j+1);
        pHi_minus = pHi(j);
        pHi_plus = pHi(j+1);
        dt = time_plus-time_minus;
        dpHi = pHi_plus-pHi_minus;
        dpHi_dt(j) = dpHi/dt;

    end

    aux1 = find(dpHi_dt == min(dpHi_dt));
    min_dpHidt = min(dpHi_dt);  % maximal rate of decline for pHi

    % Calculate Delta_pHs    
    aux2 = find(pHs == max(pHs));
    delta_pHs = max(pHs)-pH_out;
    
    %% Figure 6 (B)
    subplot(2,2,3);
    plot(d_EUS,delta_pHs, '.','MarkerSize',27)
    ylabel('(\DeltapH_S)_{max}')
    xlabel('d (\mum)')
    title('B');
    hold all
    
    %% Figure 6 (C)
    subplot(2,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    %xlim([time_aux time(end)])
    %ylim([6.99 7.22])
    xlabel('Time (sec)')
    ylabel('pH_i');
    title('C');
    hold all
    
    %% Figure 6 (D)
    subplot(2,2,4);
    plot(d_EUS,-min_dpHidt, '.','MarkerSize',27)
    ylabel('-(dpH_i/dt)_{max}')
    xlabel('d (\mum)');
    title('D');
    hold all
    
    %% Figure 6 (A)
    subplot(2,2,1);
    plot(time_plot,pHs_plot,'LineWidth',2.5,'DisplayName',name)
    %xlim([time_aux time(end)])
    %ylim([7.4993 7.515])
    xlabel('Time (sec)')
    ylabel('pH_S')
    title('A');
    hold all       
end

%legend('d = 150 \mum','d = 100 \mum','d = 50 \mum','d = 25 \mum','d = 10 \mum','d = 5 \mum','d = 1 \mum');
lgd= legend('show');
lgd.FontSize = 16;
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('Effect of the width d of the extracellular unconvected fluid');
