% © 2011 R. Occhipinti
% SPDX-License-Identifier: GPL-3.0-or-later
%
% This is the main file to run simulations of CO2 addition 
%
% This file computes:
% 1) delta_pHs = changes in steady-state pHs
% 2) delta_pHi = changes in steady-state pHi


% Setting up the model, computing the parameter structure for the rhs in
% the ODE solver
clc
clear all
%close all

ModelParametersDistr_DE_paper % Model parameters for the distributed model

[A W] = DiffusionMatrixDistr(n_in,n_out,R,R_inf, ...
                 kappa_in,kappa_out,alpha,n_buff);

Params.DiffusionMatrix = A;
Params.BoundaryVector  = W;
Params.ReactionRates   = k;
Params.BoundaryValues  = X_inf;
Params.N               = N;
Params.NumberOfBuffers = n_buff;
Params.n_out            = n_out;
Params.n_in            = n_in;
% Solving the system 
tic
tmax = 1200;
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Stats','on');
[time,X] = ode15s(@ReactionDiffusionDistrRHS,[0,tmax],X0,options,Params);
toc
% save Simulations/PCO2m_water34/Delta100um/PlayingWithPCO2m_water/simulation_CAII20_CAIV20_PmCO2_34p2_dividedby_75000_nin325.mat

% break
% pH near the membrane

n1 = (1+n_buff)*N + n_in; % one shell below membrane; n1+1 = @membrane
tf_i = 100;
tf_s = 100;
depth = 50;      %  depth of electrode inside (in microns)
depth = 1e-4*depth;  % d in centimeters
rad_in = (R/n_in)*[0:n_in];
ind_electrode = find(rad_in >= R-depth,1); % inside

set(0,'defaultaxesfontsize',16)

figure(1)
plot(time,3-log10(X(:,n1-(n_in-ind_electrode))),'LineWidth',2)
xlim([0 time(end)])
%title('Intracellular pH vs Time @ 50\mum')
xlabel('Time (s)')
ylabel('pH_i')
hold all

figure(2)
plot(time,3-log10(X(:,n1+2)),'LineWidth',2)
xlim([0 time(end)])
%ylim([7.5 7.515])
%title('Extracellular pH vs Time')
xlabel('Time (s)')
ylabel('pH_s')
hold all
% Calculations
%break
%-------------------------------------------------------
% delta_pHs
pHs = 3-log10(X(:,n1+2));
aux = find(pHs==max(pHs));
tau_p = time(aux)
delta_pHs = max(pHs)-pH_out
% figure(30)
% plot(time, pHs,'LineWidth',2)
% hold on
%plot(time_p,max(pHs),'r*')

%--------------------------------------------------------------
% delta_pHi
pHi = 3-log10(X(:,n1-(n_in-ind_electrode)));
delta_pHi = min(pHi)-pH_in
%break
%--------------------------------------------------------------

