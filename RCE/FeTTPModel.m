clc;
clear;

%some things to do still fix the for loop at the bottom. Add in the E term! Include the calculation of the boundary layer thickness from the rotation speed.

%#ok<*NUSED>
%#ok<*GVMIS>
%#ok<*INUSD>

% CV waveform
c.scan_rate = 0.2;      % V/s
c.E_start = 0.0;
c.E_end = -2.5;
c.half_cycle_time = abs(c.E_end - c.E_start)/c.scan_rate;

time = 4*c.half_cycle_time; %total time of experiment [s] to complete 2 full cycles

%Constants
c.delta = 0.1e-4; % boundary layer thickness [m]
c.mu = 8.90e-4; % viscosity of water at 25C [Pa.s]

c.T = 298.0;
c.F = 96485.333;
c.R = 8.314459848;

c.C_Fe3_i = 0.001*1000.0; %initial Fe(III) bulk concentration at t=0 in [mol/m3] units
c.C_Fe2_i = 0.000*1000.0; %initial Fe(II) bulk concentration at t=0 in [mol/m3] units
c.C_Fe1_i = 0.000*1000.0; %initial Fe(I) bulk concentration at t=0 in [mol/m3] units
c.C_Fe0_i = 0.000*1000.0; %initial Fe(0) bulk concentration at t=0 in [mol/m3] units

c.D0_Fe3 = 1.1e-010; %Diffusion coefficient of CO2 in water at 25C at infinite dilution [m/s]
c.D0_Fe2 = 6.7e-010; %Diffusion coefficient of (CO3)2- in water at 25C at infinite dilution [m/s]
c.D0_Fe1 = 4.6e-010; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]
c.D0_Fe0 = 5.7e-010; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]

c.k0_3_2 =0.00002; % rate constant Fe(III) to Fe(II) (mol/m2/s)
c.k0_2_1 = 0.00002;% rate constant Fe(II) to Fe(I) (mol/m2/s)
c.Fe1formation = 0; % Calculation of CO rate of formation at cathode surface (mol/m2/s)
c.Fe0formation = 0; % Calculation of CO rate of formation at cathode surface (mol/m2/s)

nmesh = 500; % intial mesh
m = 0; 
xmesh = linspace(0,c.delta,nmesh); 
tspan = linspace(0,time,nmesh);

sol = pdepe(m, @(x,t,u,dudx) pde(x,t,u,dudx,c), @(x) pdeic(x,c), ...
            @(xl,ul,xr,ur,t) pdebc(xl,ul,xr,ur,t,c), xmesh, tspan);
u1 = sol(:,:,1); %sol(t(i), x(j), component)
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%Fix this still change away from the forloop and also include the calculation of the global current Fe1
current_Fe3 = zeros(1,nmesh);
current_Fe2 = zeros(1,nmesh);
current_Fe1 = zeros(1,nmesh);
current_Fe0 = zeros(1,nmesh);
for i = 1:nmesh
    dx = xmesh(nmesh) - xmesh(nmesh-1);
    current_Fe3(i) = -c.F*c.D0_Fe3*(sol(i,nmesh,1)-sol(i,nmesh-1,1))/1000/dx;
    current_Fe2(i) = -c.F*c.D0_Fe2*(sol(i,nmesh,2)-sol(i,nmesh-1,2))/1000/dx;
    current_Fe1(i) = -c.F*c.D0_Fe1*(sol(i,nmesh,3)-sol(i,nmesh-1,3))/1000/dx;
    current_Fe0(i) = -c.F*c.D0_Fe0*(sol(i,nmesh,4)-sol(i,nmesh-1,4))/1000/dx;
end
global_current = -current_Fe3+current_Fe1+2*current_Fe0

figure(1);
surf(xmesh,tspan,u1/1000.0,'edgecolor','none');
xlim([0.0, c.delta]);
%title('K^{+} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.008 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(III)[mol/L]');
view(30,20);

figure(2);
surf(xmesh,tspan,u2/1000.0,'edgecolor','none');
xlim([0.0, c.delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(II)[mol/L]');
view(30,20);

figure(3);
surf(xmesh,tspan,u3/1000.0,'edgecolor','none');
xlim([0.0, c.delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(I)[mol/L]');
view(30,20);

figure(4);
surf(xmesh,tspan,u4/1000.0,'edgecolor','none');
xlim([0.0, c.delta]);
%title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
xlabel('Distance x [m]');
ylabel('Time t [s]');
zlabel('Fe(I)[mol/L]');
view(30,20);

figure(5)
plot(tspan,current_Fe3,tspan,current_Fe2,tspan,current_Fe1,tspan,current_Fe0,tspan,global_current,'LineWidth',1.5);

figure(6)
plot(E,global_current,'LineWidth',1.5);


function [c,f,s] = pde(x,t,u,dudx,const)
    c = [1; 1; 1; 1];
    f = [const.D0_Fe3*dudx(1); const.D0_Fe2*dudx(2); const.D0_Fe1*dudx(3); const.D0_Fe0*dudx(4)];
    s = [0; 0; 0; 0]; %we change this later so that the iron can react in the entire solution
end

function u0 = pdeic(x, c) 
    u0 = [c.C_Fe3_i; c.C_Fe2_i; c.C_Fe1_i; c.C_Fe0_i];
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t,c)
    E0_3_2=-0.25; %from data
    E0_2_1=-1.3;
    E0_1_0=-2.0;

    if t <= c.half_cycle_time
        E = c.E_start - c.scan_rate*t;
    elseif (t > c.half_cycle_time) && (t <= 2*c.half_cycle_time)
        E = c.E_end + c.scan_rate*(t-c.half_cycle_time);
    elseif (t > 2*c.half_cycle_time) && (t <= 3*c.half_cycle_time)
        E = c.E_start - c.scan_rate*(t-2*c.half_cycle_time);
    elseif (t > 3*c.half_cycle_time) && (t <= 4*c.half_cycle_time)
        E = c.E_end + c.scan_rate*(t-3*c.half_cycle_time);
    end

    ql = [0; 0; 0; 0];
    pl = [ul(1)-c.C_Fe3_i; ul(2)-c.C_Fe2_i; ul(3)-c.C_Fe1_i; ul(4)-c.C_Fe0_i];
    qr = [1; 1; 1; 1];
    pr = [c.k0_3_2*(ur(1)*exp(-0.5*(E-E0_3_2)*c.F/c.R/c.T)-ur(2)*exp(0.5*(E-E0_3_2)*c.F/c.R/c.T));...
            -c.k0_3_2*(ur(1)*exp(-0.5*(E-E0_3_2)*c.F/c.R/c.T)-ur(2)*exp(0.5*(E-E0_3_2)*c.F/c.R/c.T))+c.k0_3_2*(ur(2)*exp(-0.5*(E-E0_2_1)*c.F/c.R/c.T)-ur(3)*exp(0.5*(E-E0_2_1)*c.F/c.R/c.T));...
            -c.k0_3_2*(ur(2)*exp(-0.5*(E-E0_2_1)*c.F/c.R/c.T)-ur(3)*exp(0.5*(E-E0_2_1)*c.F/c.R/c.T))+c.k0_3_2*(ur(3)*exp(-0.5*(E-E0_1_0)*c.F/c.R/c.T)-ur(4)*exp(0.5*(E-E0_1_0)*c.F/c.R/c.T));...
            -c.k0_3_2*(ur(3)*exp(-0.5*(E-E0_1_0)*c.F/c.R/c.T)-ur(4)*exp(0.5*(E-E0_1_0)*c.F/c.R/c.T))];
end