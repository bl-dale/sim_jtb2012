% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% Plotting concentration profiles as a function of space for Figure 8 
% Figure 8 (A) for the Std Experiment

function [] = create_fig8_JTB(times,Xs,n,n_ins,n_out,Ns,Rs,R_infs,n_buffs,pH_outs,Pm_CO2s,sim_dir,sim_filename_base)
%strcat(sim_dir,sim_filename_base,'.mat')
%load(strcat(sim_dir,'/',sim_filename_base,'.mat'));

shades = [150/255    75/255    0;  % brown     
            0    0.5000        0;  % green
       1.0000        0         0;  % red
            0   0.7500    0.7500;  % cyan
       0.7500        0    0.7500;  % magenta
       204/255  153/255        0;  % gold
      0.2500    0.2500    0.2500]; % black
  
set(0,'DefaultAxesColorOrder',shades)
i=1;
time   = cell2mat(times(i));
X      = cell2mat(Xs(i));
n_in   = cell2mat(n_ins(i));
N      = cell2mat(Ns(i));
n_buff = cell2mat(n_buffs(i));
R      = cell2mat(Rs(i));
R_inf  = cell2mat(R_infs(i));
pH_out = cell2mat(pH_outs(i));
Pm_CO2 = cell2mat(Pm_CO2s(i))
    
r_in   = (R/n_in)*[0:n_in];
r_out  = R + ((R_inf-R)/(n_out-1))*[0:n_out-1];
r_plot = [r_in,r_out];

I1 = find(time<=0.01, 1, 'last' );
I2 = find(time<=0.25, 1, 'last' );
I3 = find(time<=5, 1, 'last' );
I4 = find(time<=35, 1, 'last' );
I5 = find(time<=200, 1, 'last' );

I = [1 I1 I2 I3 I4 I5 length(time)];
time_aux = time(I);

set(0,'defaultaxesfontsize',18)

% Zoom in

aux = find(r_in>=0.055,1);
r_in = r_in(aux:end);
r_plot = [r_in,r_out];

figure;
subplot(1,2,1);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[R-0.1,0,0.10,0.55],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,0.10,0.55],'facecolor',c,'edgecolor',c)
set(gca,'XTick',[R-0.1  R  R+0.1])
hold on
   
for ii = 1:length(I)
    
    %h1 = figure(1);
    jinit =  aux;
    jfinal = N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5)
    
    hold all
    plot([R,R],[0,1.15*max(X(I(ii),jinit:jfinal))],'k--')
    title('P_{M,CO_{2}} = 34.20 cm/s')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([R-0.1 R+0.1])
    ylim([0 1.15*max(X(I(ii),jinit:jfinal))])    
end


%% Plotting concentration profiles as a function of space for Figure 8 (B)
clc
%clear all

% Figure 8 (B) where Pm,CO2 is divided by 10,000

%load data\PmCO2_34p2_dividedby10000.mat

shades = [150/255    75/255    0;  % brown     
            0    0.5000        0;  % green
       1.0000        0         0;  % red
            0   0.7500    0.7500;  % cyan
       0.7500        0    0.7500;  % magenta
       204/255  153/255        0;  % gold
      0.2500    0.2500    0.2500]; % black
  
set(0,'DefaultAxesColorOrder',shades)
i=5;
time   = cell2mat(times(i));
X      = cell2mat(Xs(i));
n_in   = cell2mat(n_ins(i));
N      = cell2mat(Ns(i));
n_buff = cell2mat(n_buffs(i));
R      = cell2mat(Rs(i));
R_inf  = cell2mat(R_infs(i));
pH_out = cell2mat(pH_outs(i));
Pm_CO2 = cell2mat(Pm_CO2s(i))

r_in   = (R/n_in)*[0:n_in];
r_out  = R + ((R_inf-R)/(n_out-1))*[0:n_out-1];
r_plot = [r_in,r_out];

I1 = find(time<=0.01, 1, 'last' );
I2 = find(time<=0.25, 1, 'last' );
I3 = find(time<=5, 1, 'last' );
I4 = find(time<=35, 1, 'last' );
I5 = find(time<=200, 1, 'last' );

I = [1 I1 I2 I3 I4 I5 length(time)];
time_aux = time(I);
set(0,'defaultaxesfontsize',18)

% Zoom in

aux = find(r_in>=0.055,1);
r_in = r_in(aux:end);
r_plot = [r_in,r_out];

subplot(1,2,2);
c=[178/255 220/255 1];
g = [169/255 207/255 161/255];
rectangle('Position',[R-0.1,0,0.10,0.55],'facecolor',g,'edgecolor',g)
rectangle('Position',[R,0,0.10,0.55],'facecolor',c,'edgecolor',c)
set(gca,'XTick',[R-0.1  R  R+0.1])
hold on

for ii = 1:length(I)
    
    %h1 = figure(2);
    jinit =  aux;
    jfinal = N;
    plot(r_plot,X(I(ii),jinit:jfinal),'LineWidth',2.5);
    hold all
    %plot([R,R],[0,1.15*max(X(I(ii),jinit:jfinal))],'k--')
    title('P_{M,CO_{2}} = 34.20/10^4 cm/s')
    xlabel('Radius (mm)')
    ylabel('Concentration (mM)')
    xlim([R-0.1 R+0.1])
    ylim([0 1.15*max(X(I(ii),jinit:jfinal))])

end
plot([R,R],[0,1.15*max(X(I(ii),jinit:jfinal))],'k--');
legend('t = 0 s','t \approx 0.01 s','t \approx 0.25 s','t \approx 5 s','t \approx 35 s','t \approx 200 s','t = 1,200 s','Location','Southeast');
set(gcf, 'Position', get(0, 'Screensize'));   
suptitle('CO_2 profile');

