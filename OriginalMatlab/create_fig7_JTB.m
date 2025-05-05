% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% This file plots Figure 7  
% Remember: The LineWidth for pHi is 5 for thicker line and 1.5 for the
% thinner one

function [] = create_fig7_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,pH_outs,Pm_CO2s,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('Number of data file: ');
%n = length(all_data)

shades_default = [ 
         1       0.40         0;               % orange
         0.40    0.20         0;               % brown
         0       0.75        0.75;             % cyan
         0         0          1;               % blue
         0.75      0         0.75;             % magenta
         0.7500    0.7500      0;              % yellow
         0.2500    0.2500    0.2500;           % black
         1.00      0           0;              % red
         0       0.50          0];             % green 
set(0,'defaultaxesfontsize',18);
set(0,'DefaultAxesColorOrder',shades_default)
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
    Pm_CO2 = cell2mat(Pm_CO2s(i));
   
    d_EUS = (R_inf - R)*1e4;
 

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
    
    % label
    j = 34.2/Pm_CO2;
    if j==1
        name = 'P_{M,CO_{2}} = 34.2 cm/s';
    else
        name = 'P_{M,CO_{2}} = 34.2/'+string(sprintf('%.1e',j));
    end        
    
    %% Figure 7B
    for j = 1:size(time,1)-1
        time_minus = time(j);
        time_plus = time(j+1);
        pHi_minus = pHi(j);
        pHi_plus = pHi(j+1);
        dt = time_plus-time_minus;
        dpHi = pHi_plus-pHi_minus;
        dpHi_dt(j) = dpHi/dt;

    end

    % Calculate Delta_pHs
    pHs = 3-log10(X(:,n1+2));
    aux2 = find(pHs == max(pHs));
    tau_p = time(aux2);
    delta_pHs = max(pHs)-pH_out;

    subplot(2,2,3);
    semilogx(Pm_CO2,delta_pHs,'.','MarkerSize',27)
    %ylim([0 8e-3])
    ylabel('(\DeltapH_s)_{max}')
    xlabel('P_{M,CO_2} (cm/sec)')
    title('B');
    hold on
    
    %% Figure 7 (C)
    subplot(2,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([6.99 7.22])
    xlabel('Time (sec)')
    ylabel('pH_i')
    title('C');
    hold all
    
    %% Figure 7 (D)
    aux1 = find(dpHi_dt == min(dpHi_dt));
    td = time(aux1);  % time delay
    min_dpHidt = min(dpHi_dt);  % maximal rate of decline for pHi    

    subplot(2,2,4);
    semilogx(Pm_CO2,-min_dpHidt, '.','MarkerSize',27)
    %ylim([0 3.7e-3])
    ylabel('-(dpH_i/dt)_{max}')
    xlabel('P_{M,CO_2} (cm/sec)')
    title('D');
    hold on
    
    %% Figure 7 (A)
    subplot(2,2,1);
    plot(time_plot,pHs_plot,'LineWidth',2.5,'DisplayName',name)
    xlim([time_aux time(end)])
    %ylim([7.4996 7.508])
    xlabel('Time (sec)')
    ylabel('pH_S')
    title('A');
    hold all
end
lgd= legend('show');
lgd.FontSize = 8;
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('Effect of CO_2 membrane permeability');
