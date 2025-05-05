% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% Plot concentration profiles as a function of time in different interior
% shells
% This is Figure 4 of the paper

function [] = create_fig4_JTB(time,X,n_in,N,sim_dir,sim_filename_base)
strcat(sim_dir,'/',sim_filename_base,'.mat')
%load(strcat(sim_dir,'/',sim_filename_base,'.mat'))

rr = [1:20:n_in n_in n_in+1];  % Interior shells plus "plasma membrane"
set(0,'defaultaxesfontsize',20)

shades = [102/255    51/255    0;  % brown
          0    0.5000        0;  % green
     1.0000        0         0;  % red
          0   0.7500    0.7500;  % cyan
     0.7500        0    0.7500;  % magenta
    204/255  153/255        0];  % gold

set(0,'DefaultAxesColorOrder',shades)
  
time_aux = [-100];
time_plot = [time_aux'; time];

warning('off');
figure

for jj = 1:length(rr)
    
    X_CO2   = X(:,rr(jj));  % CO2 at rr as a function of time
    X_H2CO3 = X(:,N+rr(jj)); % H2CO3 at rr as a function of time
    X_HA    = X(:,2*N+rr(jj)); % HA at rr as a function of time
    X_Hplus = X(:,2*N+N+rr(jj)); % H^{+} at rr as a function of time
    X_HCO3m = X(:,2*N+N+N+rr(jj)); % HCO3^{-} at rr as a function of time
    X_Am    = X(:,2*N+N+N+N+rr(jj)); % A^{-} at rr as a function of time
    pH = 3-log10(X_Hplus);
    
    %% pH figure
    subplot(2,3,4);
    hold all
    plot(time_plot,[pH(1)'; pH],'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([7 7.23])
    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('pH')
    xlabel('Time (sec)')
    ylabel('pH');
    %ylabel('Concentration (mM)')
    
    %% CO2 figure
    subplot(2,3,1);
    hold all
    plot(time_plot,[X_CO2(1)'; X_CO2],'LineWidth',2.5)
    xlim([time_aux time(end)])
    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('CO_2')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    legend('r \approx 8 \mum','r \approx 160 \mum','r \approx 320 \mum','r \approx 480 \mum','r \approx 640 \mum','r = 650 \mum','Location','best');
    
    
    %% H2CO3 figure    
    subplot(2,3,2);

    hold all
    plot(time_plot,[X_H2CO3(1)'; X_H2CO3],'LineWidth',2.5)
    xlim([time_aux time(end)])


    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('H_2CO_3')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
    
    %% HCO3- figure
    subplot(2,3,3);    
    
    hold all
    plot(time_plot,[X_HCO3m(1)'; X_HCO3m],'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([0 3.5])
    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('HCO_3^-')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
        
    %% A- figure
    subplot(2,3,6);
    
    hold all
    plot(time_plot,[X_Am(1)'; X_Am],'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([12 15.5])
    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('A^-_1')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')
        
    %% HA figure
    
    subplot(2,3,5);
    
    hold all
    plot(time_plot,[X_HA(1)'; X_HA],'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([12 15.5])

    set(gca,'XTick',[0 200 400 600 800 1000 1200]);
    set(gca,'linewidth', 2);
    title('HA_1')
    xlabel('Time (sec)')
    ylabel('Concentration (mM)')

end
warning('on');
set(gcf, 'Position', get(0, 'Screensize'));
suptitle('Intracellular concentration-time profiles');


