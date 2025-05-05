% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% This file plots Figure 12
% 100% Immobile
 
function [] = create_fig12_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,pH_outs,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('Number of data file: ');
%n = length(all_data)


shades = [0      0      0;  % black (10um)
       0.75      0    0.75; % magenta(50um) 
     70/255    70/255  70/255;  % (150um)
     90/255   90/255 90/255;  % (250um)
     110/255   110/255 110/255;  % (350um)
     130/255  130/255  130/255; % (450um)
     160/255  160/255  160/255; % (550um)
     190/255  190/255  190/255]; % (150um)                    
figure;
set(0,'DefaultAxesColorOrder',shades)
set(0,'defaultaxesfontsize',18)
j = 7
time   = cell2mat(times(j));
X      = cell2mat(Xs(j));
n_in   = cell2mat(n_ins(j));
N      = cell2mat(Ns(j));
n_buff = cell2mat(n_buffs(j));
R      = cell2mat(Rs(j));
R_inf  = cell2mat(R_infs(j));
    
id = [10 50 150 250 350 450 550 650];
for i=1:length(id)

    n1 = (1+n_buff)*N + n_in; % one shell below membrane; n1+1 = @membrane
    tf_i = 100;
    tf_s = 100;
    depth = id(i);      % depth of electrode inside (in microns)
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

    subplot(1,2,1);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([6.89 7.22])
    xlabel('Time (sec)')
    ylabel('pH_i')
    hold all
end
title('100% immobile buffer');

%% This file plots Figure 12 (B) 
clc    

id = [10 50 150 250 350 450 550 650];
shades = [0      0      0;  % black (10um)
         204/255  153/255   0;  % gold  (50um)
         70/255    70/255  70/255;  % (150um)
         90/255   90/255 90/255;  % (250um)
         110/255   110/255 110/255;  % (350um)
         130/255  130/255  130/255; % (450um)
         160/255  160/255  160/255; % (550um)
         190/255  190/255  190/255]; % (150um)                    
set(0,'defaultaxesfontsize',18)
set(0,'DefaultAxesColorOrder',shades)

j = 5
time   = cell2mat(times(j));
X      = cell2mat(Xs(j));
n_in   = cell2mat(n_ins(j));
N      = cell2mat(Ns(j));
n_buff = cell2mat(n_buffs(j));
R      = cell2mat(Rs(j));
R_inf  = cell2mat(R_infs(j));

for i=1:length(id)
    n1 = (1+n_buff)*N + n_in; % one shell below membrane; n1+1 = @membrane
    tf_i = 100;
    tf_s = 100;
    depth = id(i);      % depth of electrode inside (in microns)
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

    subplot(1,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([6.98 7.22])
    xlabel('Time (sec)')
    ylabel('pH_i')
    hold all
end
title('90% immobile buffer');
set(gcf, 'Position', get(0, 'Screensize'));   




