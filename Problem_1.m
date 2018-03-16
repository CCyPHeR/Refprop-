clc
clear all
close all

Tc = 647.096;               % Critical Temperature in K
Pc = 22064000;              % Critical Pressure in Pa
w = 0.3443;                 % Acentric Factor for Water
R = 8.314;                  % Universal Gas Constant (m^3.Pa/K.mol)
N = 100;                    % Number of Data Points
T0 = 300;                   % Minimum Temperature of p-T diagram (K)
T1 = 0.92*Tc;               % Maximum Temperature of p-T diagram (K)
T = linspace(T0,T1,N);      % Temperature Datapoints

% Peng-Robinson Constants
a = (0.45724*(R^2)*(Tc^2))/Pc;      % (m^6.Pa/mol^2)
b = (0.07780*R*Tc)/Pc;              % (m^3/mol)
k = 0.37464 + (1.54226*w) - (0.26992*(w^2));

p = ones(N,1);                      % Initialising Trial Values of Pressure
for i = 1:N
    Diff = 2;                       % Loop Starting Criterion
    Tr = T(i)/Tc;
    alpha = (1 + k*(1-(Tr^0.5)))^2;
    while Diff>=1e-4
        A = (alpha*a*p(i))/(R^2*T(i)^2);
        B = (b*p(i))./(R*T(i));
        cubic = [1 -(1-B) (A-(2*B)-(3*(B^2))) -((A*B)-(B^2)-(B^3))];
        Z = roots(cubic);           % Finding the roots of the cubic polynomial
        ZV = max(Z(Z>=0));          % Assigning the vapour compressibility factor
        ZL = min(Z(Z>=0));          % Assigning the liquid compressibility factor
        
        %Computing fugacity coefficients
        phi_V = exp(ZV-1-log(ZV-B)-((A/(2*sqrt(2)*B))*log((ZV+((1+sqrt(2))*B))/(ZV+((1-sqrt(2))*B)))));
        phi_L = exp(ZL-1-log(ZL-B)-((A/(2*sqrt(2)*B))*log((ZL+((1+sqrt(2))*B))/(ZL+((1-sqrt(2))*B)))));
        
        % Computing fugacities
        fV = phi_V*p(i);
        fL = phi_L*p(i);
        p(i) = (p(i)*(fL/fV));      % Pressure Correction
        Diff = abs((fL/fV)-1);      % Loop Criterion
    end
    p(i) = p(i)/10^6;                                                 % Conversion from Pa to MPa
end

% NIST Data
T1 = 300:3.5030:0.92*Tc;
P = [0.0035368;0.0043345;0.0052843;0.0064094;0.0077361;0.0092933;0.011113;0.013231;0.015685;0.018519;0.021777;0.025511;
    0.029775;0.034626;0.040129;0.046350;0.053361;0.061241;0.070070;0.079936;0.090931;0.10315;0.11670;0.13169;0.14823;
    0.16644;0.18645;0.20838;0.23237;0.25856;0.28710;0.31814;0.35815;0.38837;0.42789;0.47057;0.51661;0.56617;0.61946;
    0.67667;0.73800;0.80366;0.87386;0.94882;1.0288;1.1139;1.2045;1.3008;1.4030;1.5113;1.6261;1.7475;1.8759;2.0115;
    2.1546;2.3054;2.4643;2.6315;2.8072;2.9919;3.1858;3.3893;3.6025;3.8259;4.0597;4.3043;4.5600;4.8272;5.1062;5.3974;
    5.7010;6.0176;6.3474;6.6908;7.0482;7.4201;7.8068;8.2088;8.6264;9.0602;9.5106;9.9781;10.463;10.966;11.488];

% Plotting
plot(T,p,'linewidth',1.5);
hold on
plot(T1,P,'linewidth',1.5);
grid on
xlabel('Temperature (K)'),ylabel('Pressure (MPa)'),title('p - T Diagram of Water');
legend('Peng-Robinson EOS Method','NIST Data','Location','Northwest');

% Curve Fitting according to August Equation log(P) = (A - (B/T))
P1 = log(p(26:end)');
fit_parameters = polyfit((-1./T(26:end)),P1,1);             % Computing slope of the August Equation
H_vap = fit_parameters(1)*R;                                % Computing Enthalpy of Vaporization from Slope