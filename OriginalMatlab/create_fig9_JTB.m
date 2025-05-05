% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% Figure 9 of the paper

function [] = create_fig9_JTB(times,Xs,n,n_ins,Ns,Rs,R_infs,n_buffs,CAII_flags,CAIV_flags,CAII_ins,CAIV_outs,sim_dir,sim_filename_base)
%all_data = dir(fullfile(sim_dir,'*.mat'));
%disp('AllData')
%disp(all_data)
%disp('Number of data file: ');
%n = length(all_data)

shades = [0        0.5000    0;  % green
          204/255  153/255   0;  % gold
          1.0000        0    0;  % red
          0             0    1.0000]; % blue         
  
set(0,'DefaultAxesColorOrder',shades)
set(0,'defaultaxesfontsize',18)
figure;

for i=1:n
    %load(strcat(sim_dir,'\',all_data(i).name));
    time      = cell2mat(times(i));
    X         = cell2mat(Xs(i));
    n_in      = cell2mat(n_ins(i));
    N         = cell2mat(Ns(i));
    n_buff    = cell2mat(n_buffs(i));
    R         = cell2mat(Rs(i));
    R_inf     = cell2mat(R_infs(i));
    CAII_flag = cell2mat(CAII_flags(i));
    CAII_in   = cell2mat(CAII_ins(i));
    CAIV_flag = cell2mat(CAIV_flags(i));
    CAIV_out  = cell2mat(CAIV_outs(i));

 
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
    if CAII_flag == 0
        Ai = 1;
    else
        Ai = CAII_in;
    end
    
    if CAIV_flag == 0
        Ao = 1;
    else
        Ao = CAIV_out;
    end
    
    beta = ['A_{i} = ' + string(Ai) + ' and A_{o} = ' + string(Ao)];
    
    subplot(1,2,2);
    plot(time_plot,pHi_plot,'LineWidth',2.5)  % Linewidth 1.5 or 5
    xlim([time_aux time(end)])
    %ylim([6.99 7.22])
    xlabel('Time (sec)')
    ylabel('pH_i')
    % %legend('A_{i} = A_{o} = 1 (No CA)','A_{i} = 1, A_{o} = 20','A_{i} = 20, A_{o} = 1','A_{i} = A_{o} = 20 (Std Exp)')
    hold all
    % print -depsc2 PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHi_CAs_nolegend
    % saveas(h1,'PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHi_CAs_nolegend','emf')
    % saveas(h1,'PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHi_CAs_nolegend')

    subplot(1,2,1);
    plot(time_plot,pHs_plot,'LineWidth',2.5,'DisplayName',beta)
    xlim([time_aux time(end)])
    %ylim([7.4996 7.508])
    xlabel('Time (sec)')
    ylabel('pH_S')
    %legend('A_{i} = A_{o} = 1 (No CA)','A_{i} = 1, A_{o} = 20','A_{i} = 20, A_{o} = 1','A_{i} = A_{o} = 20 (Std Exp)')
    hold all

    %print -depsc2 PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHs_CAs_nolegend
    %saveas(h2,'PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHs_CAs_nolegend','emf')
    %saveas(h2,'PCO2m_water34/Delta100um/MobileBuffers/Figure7/pHs_CAs_nolegend')
end
legend show;
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('CA activity');
