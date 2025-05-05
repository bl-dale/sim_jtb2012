% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% This is Figure 10 of the paper

function [] = create_fig10_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,Atot_ins,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('Number of data file: ');
%n = length(all_data)

shades_default = [ 0.75      0         0.75;             % magenta
                     0       0.75        0.75;             % cyan
                     1.00      0           0;              % red
                     0         0           1;              % blue
                     0       0.50          0;              % green  
                     0         0           0];             % black
set(0,'DefaultAxesColorOrder',shades_default);
set(0,'defaultaxesfontsize',18);
figure;    

for i=1:n
    %load(strcat(sim_dir,'\',all_data(i).name));
    time    = cell2mat(times(i));
    X       = cell2mat(Xs(i));
    n_in    = cell2mat(n_ins(i));
    N       = cell2mat(Ns(i));
    n_buff  = cell2mat(n_buffs(i));
    R       = cell2mat(Rs(i));
    R_inf   = cell2mat(R_infs(i));
    Atot_in = cell2mat(Atot_ins(i))
 
    n1 = (1+n_buff)*N + n_in; % one shell below membrane; n1+1 = @membrane

    depth = 50;      %  depth of electrode inside (in microns)
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
    beta = [string(round(Atot_in/27.312560103865501, 1)) + '\cdot\beta_{std}'];
    
    subplot(1,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([4.7 7.25])
    xlabel('Time (sec)')
    ylabel('pH_i')
    hold all

    % print -depsc2 'Figure10/pHi'
    % saveas(h1,'Figure10/pHi')
    % saveas(h1,'Figure10/pHi','emf')

    subplot(1,2,1);
    plot(time_plot,pHs_plot,'LineWidth',2.5,'DisplayName',beta)
    xlim([time_aux time(end)])
    %ylim([7.4996 7.508])
    
    xlabel('Time (sec)')
    ylabel('pH_s')
    hold all
    
end
%legend('1000\cdot\beta_{std}','5\cdot\beta_{std}','2\cdot\beta_{std}','\beta_{std}','\beta_{std}/2','0\cdot\beta_{std}');
legend show
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('Intracellular buffering power effect');


