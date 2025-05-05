% © 2011 D. Calvetti, R. Occhipinti, E. Somersalo
% SPDX-License-Identifier: GPL-3.0-or-later
%
% 3 Buffers

% Model parameters for setting up the experiments and the model for the
% paper with Daniela and Erkki. We will consider only a Tris oocyte at 1.5%
% CO2
%

% 1. Geometry

R      = 0.13/2;   % Radius of the cell in cm
R_inf  = 0.15/2;   % Radius of the computational domain in cm
n_in   = 80;       % Discretization inside the cell
n_out  = 100;       % Discretization outside the cell (you need to change it depending on the size of R_inf)
N      = n_in + n_out + 1;
n_buff = 3;
rad_in = (R/n_in)*[0:n_in];
rad_out  = R + ((R_inf-R)/(n_out-1))*[0:n_out];


CO2_pc = 1.5;  % CO2 percent 1.5, 5, 10

% Choose CO2 permeability (alpha) across membrane
hm = 5*10^-7;  % membrane thickness (cm)
Pm_CO2 = (1.71e-5/hm); 

% Choose % of Immobile Buffer HA1/A1
Buff_pc = 100/100;
 
% Flags to control the implementation of some features
CAII_flag = 1; % addition (1) or not (0) of carbonic anhydrase II (CAII) 
CAIV_flag = 1; % addition (1) or not (0) of carbonic anhydrase IV (CAIV) 
layer_in_mem = 0; % mobility below the membrane (vesicles) smaller

% Implement CAII and CAIV activity
A_CAII = 1;
CAII_in =A_CAII*20; % Acceleration factor for CAII activity
A_CAIV = 1;
CAIV_out = A_CAIV*20; % Acceleration factor for CAIV activity

% 2. Mobility profile
% Mobility outside the cell. Units in cm2/s
% Order: CO2, H+, H2CO3, HCO3-
% Add buffers: HA1, A1, HA2, A2, ...
% Piecewise constant background profile
kappa_out = [1.71e-5;  8.69e-5; 1.11e-5; 1.11e-5];
kappa_out = kappa_out*ones(1,n_out);

% Mobility inside the cell
kappa_in  = [1.71e-5;  8.69e-5; 1.11e-5; 1.11e-5];
kappa_in  = kappa_in*ones(1,n_in+1);

kappaBuff_in  = [0; 0]; % HA1, A1m
kappaBuff_out = [1.56e-5; 1.56e-5]; % HA1, A1m

% Mobility of buffers
if n_buff == 3
   kappaBuff_in  = [kappaBuff_in; 1.56e-5; 1.56e-5]; % Add HA2, A2m
   kappaBuff_out = [kappaBuff_out; 1.56e-5; 1.56e-5]; % Add HA2, A2m
end

kappaBuff_in = kappaBuff_in*ones(1,n_in+1);
kappaBuff_out = kappaBuff_out*ones(1,n_out);
kappa_in  = [kappa_in;kappaBuff_in];
kappa_out = [kappa_out;kappaBuff_out];

% Adding a layer of lower mobility (vesicles)
if layer_in_mem == 1
    d1 = 10;      %  distance of the layer from the surface (membrane) in microns
    d1 = 1e-4*d1;  % d in centimeters
    ind1 = find(rad_in <= R-d1);
    d2 = 50;      % depth of the layer in microns
    d2 = 1e-4*d2; % d in centimeters
    ind2 = find(rad_in < R-(d1+d2));
    ind = setdiff(ind1,ind2);
    tiny = .125;  % reduce mobility in the layer
    kappa_in(1:end,ind) = tiny*kappa_in(1:end,ind);
end

% 3. Reaction rate profiles. Order
% kb(1) : CO2 + H2O -> H2CO3  , background value
% kb(2) : H2CO3 -> CO2 + H2O
% kb(3) : H2CO3 -> HCO3- + H+
% kb(4) : HCO3- + H+ -> H2CO3
% Add buffers : k(5) HA1 -> A1- + H+
%               k(6) A1- + H+ -> HA1    etc.
% Constant background
kb = NaN(2*(1+n_buff),1);
kb(1) = 0.0302;        % Reaction rate of CO2 + H2O --> H2CO3  (s^-1)
kb(2) = 10.9631;       % Reaction velocity of H2CO3 --> CO2 + H2O  (s^-1)
K1=kb(1)/kb(2);           
pK2 = 3.618357367951740;   % Dissociation const for carbonic acid
K2   = 10^(-pK2+3);   % Units in mM
kb(3) = 1e16;         % This value is a guess chosen to satisfy the equilibrium
kb(4) = kb(3)/K2;
pK_CO2 = -log10(K1)+pK2; % overall pK for CO2/HCO3m

pKHA1_out = 7.5;
KHA1_out  = 10^(-pKHA1_out+3);
kb(5) = 1e10;
kb(6) = kb(5)/KHA1_out;

if n_buff == 3
    pKHA2_out = 7.5;
    KHA2_out  = 10^(-pKHA2_out+3);
    kb(7) = 1e10;
    kb(8) = kb(7)/KHA2_out;
end

k = kb*ones(1,N);

% Adding carbonic anhydrase II everywhere inside the oocyte
if CAII_flag == 1
    ind_in = find(rad_in<=R);
    k(1:2,ind_in) = CAII_in*k(1:2,ind_in); 
end

% Adding carbonic anhydrase IV at 5nm above the membrane
if CAIV_flag == 1
    d = .005;      %  CAIV layer in microns
    d = 1e-4*d;  % d in centimeters
    ind_out = min(n_in+find(rad_out > R+d));  %outside
    k(1:2,ind_out) = CAIV_out*k(1:2,ind_out); 
end

% 4. Transmembrane mobility

alpha      = (1/1e-20)*ones(2+2*n_buff,1);
alpha(1)   =  1/Pm_CO2; %  1/Permeability of CO2

% General Constants @ 22C (room temperature)
PB     = 760;   % mmHg
PH2O   = 35;  % mmHg
sCO2   = 0.0434; % mM/mmHg CO2 solubility @ 22C
PCO2   = CO2_pc*(PB-PH2O)/100;  % vapor pressure

% 5. Initial concentrations. Outside the cell the concentration is assumed to
% be equal to the boundary value

CO2_out    = sCO2*PCO2; %mM
pH_out     = 7.5;
Hplus_out  = 10^(-pH_out+3); % mM
H2CO3_out  = K1*CO2_out;
HCO3m_out  = (K2*H2CO3_out)/Hplus_out;

Atot1_out   = 5; % mM
HA1_out     = Hplus_out*Atot1_out/(KHA1_out + Hplus_out);
A1m_out     = Atot1_out - HA1_out;
    
if n_buff == 3
    Atot2_out   = 0; % mM
    HA2_out     = Hplus_out*Atot2_out/(KHA2_out + Hplus_out);
    A2m_out     = Atot2_out - HA2_out;
end

CO2_in = 0;
pH_in = 7.2;
        
Hplus_in  = 10^(-pH_in+3);
H2CO3_in  = K1*CO2_in;
HCO3m_in  = K2*H2CO3_in/Hplus_in;

pKHA1_in = 7.10;  % (Corresponds to a pH_fin = 7.00)
KHA1_in = 10^(-pKHA1_in+3);

kb_HA1_in_plus = 1e10;
kb_HA1_in_minus = kb_HA1_in_plus/KHA1_in;
k(5,1:n_in+1) = kb_HA1_in_plus;
k(6,1:n_in+1) = kb_HA1_in_minus;

Atot_in = 27.312560103865501; % mM  (Buffer power of 15.65 mM/pH)

if n_buff == 2
    Atot1_in = Atot_in; %mM  
elseif n_buff == 3
    Atot1_in = Buff_pc*Atot_in; % Immobile Buffer
end

HA1_in     = Hplus_in*Atot1_in/(KHA1_in + Hplus_in);
A1m_in     = Atot1_in - HA1_in;
    
pKHA2_in = 7.10;  % (Corresponds to a pH_fin = 7.0)
KHA2_in = 10^(-pKHA2_in+3);

if n_buff == 3
    
    kb_HA2_in_plus = 1e10;
    kb_HA2_in_minus = kb_HA2_in_plus/KHA2_in;
    k(7,1:n_in+1) = kb_HA2_in_plus;
    k(8,1:n_in+1) = kb_HA2_in_minus;
    
    Atot2_in = (1-Buff_pc)*Atot_in; %mM  Mobile Buffer
    HA2_in     = Hplus_in*Atot2_in/(KHA2_in + Hplus_in);
    A2m_in     = Atot2_in - HA2_in;
    
end

% Setting up an initial value vector

u0_in  = CO2_in*ones(n_in+1,1);
u0_out = CO2_out*ones(n_out,1);
u1_in  = H2CO3_in*ones(n_in+1,1);
u1_out = H2CO3_out*ones(n_out,1);
v0_in  = Hplus_in*ones(n_in+1,1);
v0_out = Hplus_out*ones(n_out,1);
v1_in  = HCO3m_in*ones(n_in+1,1);
v1_out = HCO3m_out*ones(n_out,1);

u2_in  = HA1_in*ones(n_in+1,1);
u2_out = HA1_out*ones(n_out,1);
v2_in  = A1m_in*ones(n_in+1,1);
v2_out = A1m_out*ones(n_out,1);

if n_buff == 3
    u3_in  = HA2_in*ones(n_in+1,1);
    u3_out = HA2_out*ones(n_out,1);
    v3_in  = A2m_in*ones(n_in+1,1);
    v3_out = A2m_out*ones(n_out,1);
end

u0     = [u0_in;u0_out]; % CO2
u1     = [u1_in;u1_out]; % H2CO3
v0     = [v0_in;v0_out]; % Hplus
v1     = [v1_in;v1_out]; % HCO3m

u2     = [u2_in;u2_out]; % HA1
v2     = [v2_in;v2_out]; % A1m

if n_buff == 3
    u3     = [u3_in;u3_out]; % HA2
    v3     = [v3_in;v3_out]; % A2m
end

u = [u0;u1];
v = [v0;v1];
u = [u;u2];
v = [v;v2];

if n_buff == 3
    u = [u;u3];
    v = [v;v3];
end

X0 = [u;v];

% Boundary values
u_inf = [u0_out(1);u1_out(1)];
v_inf = [v0_out(1);v1_out(1)];
u_inf = [u_inf;u2_out(1)];
v_inf = [v_inf;v2_out(1)];
if n_buff == 3
    u_inf = [u_inf;u3_out(1)];
    v_inf = [v_inf;v3_out(1)];
end

X_inf = [u_inf;v_inf];

