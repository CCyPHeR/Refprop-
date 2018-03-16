%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
%
% Example file for the Matlab R2006a+ version of the NIST Refprop interface
%
% written by J. Lux, 2007-02-22, johannes.lux@dlr.de
% 
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
% refpropm Thermophysical properties of pure substances and mixtures.
%   Calling sequence for pure substances:
%      result=refpropm(prop_req, spec1, value1, spec2, value2, substance1)
%   and for mixtures
%       result=refpropm(prop_req, spec1, value1, spec2, value2, substance1, substance2, ...,x)
%   where
%       prop_req    is a character string showing what properties that are requested
%                   Each property is represented by one character:
%                           P   Pressure [kPa]
%                           T   Temperature [K]
%                           D   Density [kg/m3]
%                           H   Enthalpy [J/kg]
%                           S   Entropy [J/(kg/K)]
%                           U   Internal energy [J/kg]
%                           C   Cp [J/(kg K)]
%                           O   Cv [J/(kg K)]
%                           K   Ratio of specific heats (Cp/Cv) [-]
%                           A   Speed of sound [m/s]
%                           X   liquid phase and gas phase composition (mass fractions)
%                           V   Dynamic viscosity [Pa*s]
%                           L   Thermal conductivity [W/(m K)]
%                           Q   Quality (vapor fraction) (kg/kg)
%                           I   Surface tension [N/m]
%
%       spec1           is a character giving what we want to specify (T, P, H or D)
%       value1          is the corresponding value
%       spec2           is a character giving the second specification (P, D, H, S, U or Q)
%       value2          is the value of the second specification
%
%       substance1  is a string with the name of the first (maybe only) substance
%       substance2,..., substance N  are the name of the other substances in the mixture.
%                       Up to 20 substances can be handled
%                       Valid substance names are equal to the file names in the
%                       C:\Program Files\REFPROP\fluids\' directory (with .FLD excluded).
%                       is a vector with mass fractions of the substances in the mixture.
%
%   Examples:
%   1) P=refpropm('P','T',373.15,'Q',0,'water') gives
%         the vapor pressure of water at 373.15 K in [kPa]
%   2) [S Cp]=refpropm('SC','T',373.15,'Q',1,'water') gives
%         Entropy and Cp of saturated steam at 373.15 K
%   3) viscmix=refpropm('V','T',323.15,'P',1e5,'water','ammonia',[0.9 0.1]) gives
%         the viscosity of a 10% ammonia in water at 100 kPa and 323.15 K.
%   4) [x y]=refpropm('X','P',5e5,'Q',0.4,'R134a','R32',[0.8, 0.2]) gives
%         temperature as well as gas and liquid compositions for a mixture
%         of two refrigerants at a certain pressure and quality.
%         Note that two output variables are needed when 'X' is requested.
%
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
clear all;
clc;
% calculation of the density of methane at constant temperature
temp=200;           % [K]
fluid='methane';

for p=1:10      % 1...10 bar
    p=p*1e2;    % input pressure unit is [kPa]
    rho_ch4=refpropm('D','T',temp,'P',p,fluid);     % fluid property calculation
    disp(['rho=',num2str(rho_ch4),' kg/m³'])
end
