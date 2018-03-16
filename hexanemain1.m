%%Peng Robinson - hexane
clc;
clear all;

accfact = 0.299;
Tc = 507.82; %% in kelvin
Pc = 3.034; %% in Mpa
Dc = 233.18000000; %% kg/m3

%MW = 86.18*10^(-3);
R = 0.008314; %% J/mol-K
K = 0.37464 + 1.54226*accfact - 0.26992*(accfact^(2));
a = (0.45724*R^2*Tc^2)/Pc;
%a = (0.42748*R^2*Tc^2)/Pc;
b = 0.0778*R*Tc/Pc;

filename2='G:\Masters-TUDelft\Year1\Sem2\AppliedThermo\fluidhexane13p.xlsx';
%density=xlsread(filename2,'C2:C12');
T = xlsread(filename2,'A2:A14');
%Vm = xlsread(filename2,'D2:D14');
% i=2;
% P(i) = 0.01;

for i =1:13
%Vm(i) = MW/density(i);
Tr = T(i)/Tc;
alpha = (1+K*(1-Tr^(0.5)))^2;

if i==1
    pg=0.01;
else
    pg = P(i-1);
end

diff=1.1;

while abs(diff-1)>10^(-4)

pg = pg*diff;
A = a*alpha*pg/(R^2*T(i)^2);

B = b*pg/(R*T(i));

coeff = [1, -(1-B), (A-2*B-3*B*B), -A*B+B^2+B^3];
Z = roots(coeff);

Zl = min(Z);
Zv = max(Z);

fl = exp(Zl-1-log(Zl-B)-(A/(2^(1.5)*B))*log((Zl+(1+2^(0.5))*B)/(Zl+(1-2^(0.5))*B)));
fv = exp(Zv-1-log(Zv-B)-(A/(2^(1.5)*B))*log((Zv+(1+2^(0.5))*B)/(Zv+(1-2^(0.5))*B)));

diff = (fl/fv);

% P(i) = (R*T(i)/(Vm(i)-b)) - ((a*alpha)/(Vm(i)^2+2*b*Vm(i)-b^2));
% 
% paug(i) = exp(14.0568-2825.42/(T(i)-42.7089));

end
P(i) = pg;
end

% pnist = xlsread(filename2,'B2:B14');
% Tinv = 1./T;
% 
% hold on;
% plot(T,pnist,'-r','LineWidth',2);
% plot(T,P,'bd','MarkerSize',10);
% 
% figure(1)