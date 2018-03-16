clear
clc

N = 100;                % Number of Data Points
x1 = linspace(0,1,N);   % Mole Fraction of Ethanol in Liquid Phase
x2 = 1-x1;              % Mole Fraction of Hexane in Liquid Phase

A12 = 0.0952;           % Wilson Model Parameters
A21 = 0.2713;
P = 101.3;              % Total Pressure in kPa

%Calculating boiling points of pure components
Tb1 = 55.7172 + (3423.53/(16.1950 - log(P)))-273.15; % Boiling Point of Ethanol in K
Tb2 = 42.7089 + (2825.42/(14.0568 - log(P)))-273.15; % Boiling Point of Hexane in K

% Calculating Activity Coefficients with Wilson Model
gamma1 = exp((-log(x1+(x2.*A12))) + (x2.*((A12./(x1+(x2.*A12))) - (A21./(x2+(x1.*A21))))));
gamma2 = exp((-log(x2+(x1.*A21))) + (x1.*((A21./(x2+(x1.*A21))) - (A12./(x1+(x2.*A12))))));

% Solving the system of Equations for T
Equation = @(T) P - (gamma2.*exp(14.0568-(2825.42./((T+273.15)-42.7089))))- x1.*((gamma1 ...
       .*exp(16.1952-(3423.53./((T+273.15)-55.7172))))-(gamma2.*exp(14.0568-(2825.42./((T+273.15)-42.7089)))));
    
T(1:N) = fsolve(Equation,1:N);

%Calculating Saturation Pressures
p1 = exp(16.1952-(3423.53./((T+273.15)-55.7172))); % T in Celsius and p1 in kPa
p2 = exp(14.0568-(2825.42./((T+273.15)-42.7089)));

% Calculating mole fractions in vapour phase
y1 = (x1.*gamma1.*p1)/P;    % Mole Fraction of Ethanol in Vapor phase
y2 = 1 - y1;                % Mole Fraction of Hexane in Vapor Phase

% Plotting T-xy Diagram
[x,y] = intersections(x1(2:end-1),T(2:end-1),y1(2:end-1),T(2:end-1),1);% Finding the coordinate of azeotrope formation
plot(x1,T,'linewidth',1.5);
hold on
plot(y1,T,'linewidth',1.5);
plot(x,y,'c.','markersize',12);
xlabel('Mole Fraction of Ethanol'),ylabel('Temperature [^oC]');
title('T-xy Diagram for Ethanol/Hexane at 101.3 kPa');
legend('Mole Fraction in Liquid Phase','Mole Fraction in Vapour Phase','Azeotrope Formation Point','Location','Northwest');

% Available Data from Dortmund Data Bank

X1 = [0;0.01;0.02;0.06;0.08;0.152;0.245;0.333;0.452;0.588;0.670;0.725;0.765;0.898;0.955;0.99;0.994;1];
Y1 = [0;0.095;0.193;0.365;0.42;0.532;0.605;0.63;0.64;0.65;0.66;0.67;0.675;0.71;0.745;0.84;0.935;1];
X2 = 1-X1;              % Mole Fraction of Ethanol in Liquid Phase
Y2 = 1-Y1;              % Mole Fraction of Ethanol in Vapor Phase
T1 = [351.45,349.15,346.35,340.55,339.05,334.95,332.55,331.85,331.5,331.25,331.15,331.4,331.6,332.3,333.35, ...
    336.65,339.85,341.85];
T1 = T1-273.15;