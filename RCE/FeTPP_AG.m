%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%        1D model for calculation of FeTPP CVs                              %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%                Morales-Guio Group. UCLA                                   %%%%%%%%%%%%%
%%%%%%%%%%%%%                                                                           %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

clc;
clear;

global delta nmesh 
global scan_rate E_start E_end potential_range half_cycle_time

% Generate CV waveform
scan_rate = 0.2      % V/s
E_start = 0.0
E_end = -2.5
potential_range = abs(E_end - E_start)
half_cycle_time = potential_range/scan_rate

time = 4*half_cycle_time; %total time of experiment [s] to complete 2 full cycles

%------------------------------
%Deffinition of Constants
%------------------------------
thickness =0.1E-4;% boundary layer thickness [m]

% cefunknown;
delta = thickness; % boundary layer thickness [m]
mu = 8.90e-004; % viscosity of water at 25C [Pa.s]

global T R F
T = 298.0; %Temperature [K]

%------------------------------
%Physical constants
F = 96485.333; %Faraday's constant [C/mol of e-]
R = 8.314459848; %Gas constant [J/K/mol]
%------------------------------ 
%electrochemical reactions on the electrode
% 1 electron reactions
%   Fe(III) + 1e- <-> Fe(II)
%   Fe(II) + 1e- <-> Fe(I)
%   Fe(I) + 1e- <-> Fe(0)



global C_Fe3_i C_Fe2_i C_Fe1_i C_Fe0_i 
C_Fe3_i = 0.001*1000.0; %initial Fe(III) bulk concentration at t=0 in [mol/m3] units
C_Fe2_i = 0.000*1000.0; %initial Fe(II) bulk concentration at t=0 in [mol/m3] units
C_Fe1_i = 0.000*1000.0; %initial Fe(I) bulk concentration at t=0 in [mol/m3] units
C_Fe0_i = 0.000*1000.0; %initial Fe(0) bulk concentration at t=0 in [mol/m3] units
%-----------------------------
%Diffusion coefficients 
global D0_Fe3 D0_Fe2 D0_Fe1 D0_Fe0
D0_Fe3 = 1.1e-010; %Diffusion coefficient of CO2 in water at 25C at infinite dilution [m/s]
D0_Fe2 = 6.7e-010; %Diffusion coefficient of (CO3)2- in water at 25C at infinite dilution [m/s]
D0_Fe1 = 4.6e-010; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]
D0_Fe0 = 5.7e-010; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]

%-----------------------------
%Diffusion coefficients corrected for varying electrolyte concentration
%using Stokes-Einsteinn's equation

% %Ionic mobilities corrected for varying electrolyte concentration
% u_CO32m = D_CO32m %mobility of (CO3)2- in HCO3- solution at 25C [m/s]
% u_HCO3m = D_HCO3m; %mobility of HCO3- in HCO3- solution at 25C [m/s]
% u_OHm = D_OHm; %mobility of OH- in HCO3- solution at 25C [m/s]
% u_HCOOm = D_HCOOm; %mobility of formate in HCO3- solution at 25C [m/s]
% u_acetate = D_acetate; %mobility of acetate in HCO3- solution at 25C [m/s]
% u_K = D_K %mobility of potassium ion in HCO3- solution at 25C [m/s]
%-----------------------------
global k0_3_2 k0_2_1 Fe1formation Fe0formation 

% rate constant Fe(III) to Fe(II) (mol/m2/s)
k0_3_2 =0.00002;

% rate constant Fe(II) to Fe(I) (mol/m2/s)
k0_2_1 = 0.00002;
% Calculation of CO rate of formation at cathode surface (mol/m2/s)
Fe1formation = 0;
% Calculation of CO rate of formation at cathode surface (mol/m2/s)
Fe0formation = 0;

vpot_i = -0.6;



%-----------------------------
%SOLUTION BLOCK
%-----------------------------

nmesh = 501; % intial mesh
xplot = nmesh; % meshes for plotting the result on x.
tplot = nmesh; % meshes for plotting the result on t.
% solution block.
m = 0; 
x = linspace(0,delta,nmesh); 
t = linspace(0,time,nmesh);
t_grid = linspace(0,time,nmesh);


for l=1:nmesh
    if t_grid(l) <= half_cycle_time
     E(l) = E_start - scan_rate*t_grid(l);
    elseif (t_grid(l) > half_cycle_time) && (t_grid(l) <= 2*half_cycle_time)
      E(l) = E(l-1) + scan_rate*(time/(nmesh-1));
    elseif (t_grid(l) > 2*half_cycle_time) && (t_grid(l) <= 3*half_cycle_time)
      E(l) = E(l-1) - scan_rate*(time/(nmesh-1));
    elseif (t_grid(l) > 3*half_cycle_time) && (t_grid(l) <= 4*half_cycle_time)
      E(l) = E(l-1) + scan_rate*(time/(nmesh-1));
    end
end

E;
sol = pdepe(m,@pdex5pde,@pdex5ic,@pdex5bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);

%0.05916*(14+log10(abs(sol(:,:,4)/1000))-pH_i)

% Compute individual currents
current_Fe3 = zeros(1,nmesh);
current_Fe2 = zeros(1,nmesh);
current_Fe1 = zeros(1,nmesh);
current_Fe0 = zeros(1,nmesh);
for i = 1:nmesh
    dx = x(nmesh) - x(nmesh-1);
    current_Fe3(i) = -F*D0_Fe3*(sol(i,nmesh,1)-sol(i,nmesh-1,1))/1000/dx;
    current_Fe2(i) = -F*D0_Fe2*(sol(i,nmesh,2)-sol(i,nmesh-1,2))/1000/dx;
    current_Fe1(i) = -F*D0_Fe1*(sol(i,nmesh,3)-sol(i,nmesh-1,3))/1000/dx;
    current_Fe0(i) = -F*D0_Fe1*(sol(i,nmesh,4)-sol(i,nmesh-1,4))/1000/dx;
end
x(nmesh)
x(nmesh-1)
global_current = -current_Fe3+current_Fe1+2*current_Fe0

% % Generate corresponding E(t)
% scan_rate = 0.05;
% E_start = 0.0;
% E_end = -2.5;
% potential_range = abs(E_end - E_start);
% half_cycle_time = potential_range / scan_rate;
% E = zeros(nmesh);
% for l=1:nmesh
%     if t_grid <= half_cycle_time
%         E(l) = E_start - scan_rate*t_grid;
%     elseif t <= 2*half_cycle_time
%         E(l) = E_end + scan_rate*(t_grid-half_cycle_time);
%     elseif t <= 3*half_cycle_time
%         E(l) = E_start - scan_rate*(t_grid-2*half_cycle_time);
%     elseif t <= 4*half_cycle_time
%         E(l) = E_end + scan_rate*(t_grid-3*half_cycle_time);
%     end
% end

% Plot CV
% figure;
% plot(E, global_current, 'k', E, current_Fe3, 'r--', E, current_Fe2, 'b--', E, current_Fe1, 'g--');
% xlabel('Potential (V vs. Ag/AgCl)');
% ylabel('Current (A/m^2)');
% title('Simulated CV of FeTPP');
% legend('Total Current', 'Fe^{3+} Flux', 'Fe^{2+} Flux', 'Fe^{1+} Flux');
% grid on;
% set(gca, 'XDir', 'reverse');


Fe3conc = sol(nmesh,:,1)/1000.0;
Fe2conc = sol(nmesh,:,2)/1000.0;
Fe1conc = sol(nmesh,:,3)/1000.0;
Fe0conc = sol(nmesh,:,4)/1000.0;


figure(1);
 surf(x,t,u1/1000.0,'edgecolor','none');
 set(gca,'xlim',[0.0 delta]);
 set(gca,'fontsize',11.5);
 %title('K^{+} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.008 cm');
 xlabel('Distance x [m]');
 ylabel('Time t [s]');
 zlabel('Fe(III)[mol/L]');
 view(30,20);



figure(2);
surf(x,t,u2/1000.0,'edgecolor','none');
set(gca,'xlim',[0.0 delta]);
set(gca,'fontsize',11.5);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(II)[mol/L]');
view(30,20);

figure(3);
surf(x,t,u3/1000.0,'edgecolor','none');
set(gca,'xlim',[0.0 delta]);
set(gca,'fontsize',11.5);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(I)[mol/L]');
view(30,20);

figure(4);
surf(x,t,u4/1000.0,'edgecolor','none');
set(gca,'xlim',[0.0 delta]);
set(gca,'fontsize',11.5);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(I)[mol/L]');
view(30,20);


% figure
% %title({'Steady state CO_{2}, HCOO^{-}, OH^{-} conc. in KHCO_{3} = 0.1 M,\delta = 0.01 cm';' '});
% hold on
% 
% yyaxis right
% plot(x,Fe3conc,'LineWidth',1.5);
% ylabel('Fe^{III} concentration [M]');
% set(gca,'xlim',[0.0 delta]);
% set(gca,'fontsize',11.5);
% xlabel('Distance x [m]');
% 
% yyaxis left
% plot(x,Fe3conc,x,Fe2conc,x,Fe1conc,'LineWidth',1.5);
% ylabel('Concentration [M]');
% 
% hold off
% legend('Fe3conc','Fe2conc','Fe3conc','Location','west');

figure(5)
plot(t,current_Fe3,t,current_Fe2,t,current_Fe1,t,current_Fe0,t,global_current,'LineWidth',1.5);

figure(6)
plot(E,global_current,'LineWidth',1.5);

%title({'Steady state pH in KHCO_{3} = 0.1 M, \delta = 0.01 cm';' '});


%final concentrations on the electrode surface
u1final = u1(nmesh,nmesh)/1000.0;
u2final = u2(nmesh,nmesh)/1000.0;
u3final = u3(nmesh,nmesh)/1000.0;
u4final = u4(nmesh,nmesh)/1000.0;


%------------------------------
%This function below extrapolates viscosity from experimental
%the viscosity.txt file must be saved in the same folder as this file
%------------------------------



function [c,f,s] = pdex5pde(x,t,u,DuDx)
%Diffusion coefficients corrected for varying electrolyte concentration
%using Stokes-Einsteinn's equation

global D0_Fe3 D0_Fe2 D0_Fe1 D0_Fe0 
global e
e = 1.602176620898E-19;
%-----------------------------
c = [1; 1; 1; 1];
f = [D0_Fe3*DuDx(1); D0_Fe2*DuDx(2); D0_Fe1*DuDx(3); D0_Fe0*DuDx(4)];
s1 = 0;
s2 = 0;
s3 = 0;
s4 = 0;
s = [s1; s2; s3; s4];
end
%u(16)-u(2)-2.0*u(3)-u(4)-u(9)-u(14)
% --------------------------------------------------------------------------

function u0 = pdex5ic(x) 
% Initial composition of the bulk electrolyte at t=0
global C_Fe3_i C_Fe2_i C_Fe1_i C_Fe0_i 

%------------------------------
u0 = [C_Fe3_i; C_Fe2_i; C_Fe1_i; C_Fe0_i];
end

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex5bc(xl,ul,xr,ur,t)
% Initial composition of the bulk electrolyte at t=0
global C_Fe3_i C_Fe2_i C_Fe1_i C_Fe0_i R T F
global k0_3_2 k0_2_1 Fe1formation Fe0formation 
global scan_rate E_start E_end half_cycle_time


%------------------------------
pl = [ul(1)-C_Fe3_i; ul(2)-C_Fe2_i; ul(3)-C_Fe1_i; ul(4)-C_Fe0_i];
ql = [0; 0; 0; 0];

E0_3_2=-0.25;
E0_2_1=-1.3;
E0_1_0=-2.0;

% % Generate CV waveform
% scan_rate = 0.05;      % V/s
% E_start = 0.0;
% E_end = -2.5;
% potential_range = abs(E_end - E_start);
% half_cycle_time = potential_range/scan_rate;

if t <= half_cycle_time
     E = E_start - scan_rate*t;
elseif (t > half_cycle_time) && (t <= 2*half_cycle_time)
      E = E_end + scan_rate*(t-half_cycle_time);
elseif (t > 2*half_cycle_time) && (t <= 3*half_cycle_time)
      E = E_start - scan_rate*(t-2*half_cycle_time);
elseif (t > 3*half_cycle_time) && (t <= 4*half_cycle_time)
      E = E_end + scan_rate*(t-3*half_cycle_time);
end

% t
% E
pr = [k0_3_2*(ur(1)*exp(-0.5*(E-E0_3_2)*F/R/T)-ur(2)*exp(0.5*(E-E0_3_2)*F/R/T));...
    -k0_3_2*(ur(1)*exp(-0.5*(E-E0_3_2)*F/R/T)-ur(2)*exp(0.5*(E-E0_3_2)*F/R/T))+k0_3_2*(ur(2)*exp(-0.5*(E-E0_2_1)*F/R/T)-ur(3)*exp(0.5*(E-E0_2_1)*F/R/T));...
    -k0_3_2*(ur(2)*exp(-0.5*(E-E0_2_1)*F/R/T)-ur(3)*exp(0.5*(E-E0_2_1)*F/R/T))+k0_3_2*(ur(3)*exp(-0.5*(E-E0_1_0)*F/R/T)-ur(4)*exp(0.5*(E-E0_1_0)*F/R/T));...
    -k0_3_2*(ur(3)*exp(-0.5*(E-E0_1_0)*F/R/T)-ur(4)*exp(0.5*(E-E0_1_0)*F/R/T))];
qr = [1; 1; 1; 1];
end
