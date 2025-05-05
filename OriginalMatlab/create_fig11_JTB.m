% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% This file plots Figure 11 (A) and (B)

function [] = create_fig11_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,Buff_pcs,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('Number of data file: ');
%n = length(all_data)

shades = [0         0        1;  % blue
          0      0.50        0;  % green 
          1      0.40        0;  % orange
       0.40      0.20        0;  % brown
    204/255   153/255        0;  % gold
          0   153/255  204/255;  % cyan
       0.75         0     0.75]; % magenta                            
figure;
set(0,'DefaultAxesColorOrder',shades)
set(0,'defaultaxesfontsize',18)

for i=1:n
    %load(strcat(sim_dir,'\',all_data(i).name));
 
    time    = cell2mat(times(i));
    X       = cell2mat(Xs(i));
    n_in    = cell2mat(n_ins(i));
    N       = cell2mat(Ns(i));
    n_buff  = cell2mat(n_buffs(i));
    R       = cell2mat(Rs(i));
    R_inf   = cell2mat(R_infs(i));
    Buff_pc = cell2mat(Buff_pcs(i));
    
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
    name = string(Buff_pc)+ '%';
    
    subplot(1,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)
    xlim([time_aux time(end)])
    %ylim([6.92 7.22])
    % title('\delta = 100 \mum')
    xlabel('Time (sec)')
    ylabel('pH_i')
    % legend('0%','50%','70%','80%','90%','95%','100%')
    hold all
    % print -depsc2 PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHi_ThreeBuffers
    % saveas(h1,'PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHi_ThreeBuffers','emf')
    % saveas(h1,'PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHi_ThreeBuffers')

    subplot(1,2,1);
    plot(time_plot,pHs_plot,'LineWidth',2.5,'DisplayName',name)
    xlim([time_aux time(end)])
    %ylim([7.4996 7.508])
    %  title('\delta = 100 \mum')
    xlabel('Time (sec)')
    ylabel('pH_S')
    
    % vv = get(hh,'title');
    % set(vv,'string','Fraction of total buffer\newline that is immobile ([TA_2]_i)');
    hold all

    % print -depsc2 PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHs_ThreeBuffers
    % saveas(h2,'PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHs_ThreeBuffers','emf')
    % saveas(h2,'PCO2m_water34/Delta100um/MobileBuffers/Figure11/pHs_ThreeBuffers')
    % 
end
%hh = legend('0%','50%','70%','80%','90%','95%','100%');
leg = legend('show');
title(leg,{'Fraction of total buffer','that is immobile ([TA_2]_i)'});
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('Mobile versus immobile intracellular buffers');
