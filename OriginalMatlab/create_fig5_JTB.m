% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% Plotting concentration profiles as a function of distance from the
% center
% This is Figure 5 of the paper

function [] = create_fig5_JTB(time,X,n_in,n_out,N,R,R_inf,sim_dir,sim_filename_base)
%strcat(sim_dir,'/',sim_filename_base,'.mat')
%load(strcat(sim_dir,'/',sim_filename_base,'.mat'))

r_in   = (R/n_in)*[0:n_in];
r_out  = R + ((R_inf-R)/(n_out-1))*[0:n_out-1];
r_plot = [r_in,r_out];

I = [1:610:length(time) length(time)];
time_aux = time(I);

set(0,'defaultaxesfontsize',20);
warning('off');
figure
shades = [102/255    51/255  0;  % brown
          0    0.5000        0;  % green
     1.0000        0         0;  % red
          0   0.7500    0.7500;  % cyan
     0.7500        0    0.7500;  % magenta
    204/255  153/255        0];  % gold
  
set(0,'DefaultAxesColorOrder',shades)

%% color for background
subplot(2,3,1);
c = [178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[0, 0, R, max(1.3*max(X(:,1:N)))],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,R_inf-R,max(1.3*max(X(:,1:N)))],'facecolor',c,'edgecolor',c)
hold on

subplot(2,3,2);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[0, 0, R, max(1.3*max(X(:,1+N:2*N)))],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,R_inf-R,max(1.3*max(X(:,1+N:2*N)))],'facecolor',c,'edgecolor',c)
hold on

subplot(2,3,5);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[0, 0, R, max(1.3*max(X(:,1+2*N:3*N)))],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,R_inf-R,max(1.3*max(X(:,1+2*N:3*N)))],'facecolor',c,'edgecolor',c)
hold on

subplot(2,3,4);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
pH_1 = min(0.95*min(3-log10(X(:,1+3*N:4*N))));
pH_2 = max(1.025*max(3-log10(X(:,1+3*N:4*N))));
rectangle('Position',[0, pH_1, R, pH_2-pH_1],'facecolor',g,'edgecolor',g)
rectangle('Position',[R, pH_1, R_inf-R, pH_2-pH_1],'facecolor',c,'edgecolor',c)
hold on

subplot(2,3,3); 
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[0, 0, R, max(1.3*max(X(:,1+4*N:5*N)))],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,R_inf-R,max(1.3*max(X(:,1+4*N:5*N)))],'facecolor',c,'edgecolor',c)
hold on

subplot(2,3,6);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[0, 0, R, max(1.3*max(X(:,1+5*N:6*N)))],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,R_inf-R,max(1.3*max(X(:,1+5*N:6*N)))],'facecolor',c,'edgecolor',c)
hold on


%% plot data
for ii = 1:length(I)
    
    % CO2 figure
    subplot(2,3,1);    
    jinit =  1;
    jfinal = N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    plot([R,R],[0,1.25*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
    title('CO_2')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([0 max(r_plot)])
    ylim([0 1.3*max(X(I(ii),jinit:jfinal))])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
    
    
    % H2CO3 figure    
    subplot(2,3,2);    
    jinit  = jinit + N;
    jfinal = jfinal + N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    plot([R,R],[0,1.25*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
    title('H_2CO_3')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([0 max(r_plot)])
    ylim([0 1.3*max(X(I(ii),jinit:jfinal))])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
    
    % HA figure    
    subplot(2,3,5);        
    jinit  = jinit + N;
    jfinal = jfinal + N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    plot([R,R],[0,1.3*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
    title('HA_1')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([0 max(r_plot)])
    ylim([0 1.3*max(X(I(ii),jinit:jfinal))])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
        
    % pH figure
    subplot(2,3,4);
    
    jinit  = jinit + N;
    jfinal = jfinal + N;
    plot(r_plot,3 -log10(X(I(ii),jinit:jfinal)),'LineWidth',2.5);
    hold all
    plot([R,R],[pH_1,1.025*max(3 -log10(X(I(ii),jinit:jfinal)))],'k--','LineWidth',2)
    title('pH')
    xlabel('Radius (mm)')
    ylim([pH_1 pH_2])
    xlim([0 max(r_plot)])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
    
    % HCO3- figure
    subplot(2,3,3); 
   
    jinit  = jinit + N;
    jfinal = jfinal + N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    %plot([R,R],[0,1.22*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
    %legend('t = 0 s','t \approx 0.00005 s','t \approx 0.007 s','t \approx 0.26 s','t \approx 10.77 s','t \approx 225.08 s','t = 1,200 s','Location','best')
    title('HCO_3^{-}')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([0 max(r_plot)])
    ylim([0 1.3*max(X(I(ii),jinit:jfinal))])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
    
    % A- figure
    subplot(2,3,6);    
    jinit  = jinit + N;
    jfinal = jfinal + N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    plot([R,R],[0,1.3*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
    title('A^{-}_1')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([0 max(r_plot)])
    ylim([0 1.3*max(X(I(ii),jinit:jfinal))])
    %set(gca,'XTick',[0 0.02 0.04 0.06 0.075])
    
    %keyboard
end

%% legend
for ii = 1:length(I)
    legend_time(ii) = round(time(I(ii)),4);
end

subplot(2,3,3);
plot([R,R],[0,1.22*max(X(I(ii),jinit:jfinal))],'k--','LineWidth',2)
legend('t = ' + string(legend_time) +' sec','Location','best');
%legend('t = 0 s','t \approx 0.00005 s','t \approx 0.26 s','t \approx 10.77 s','t \approx 225.08 s','t = 1,200 s','Location','best')

warning('on');
set(gcf, 'Position', get(0, 'Screensize'));
suptitle('Extra- and intracellular concentration-distance profiles');
