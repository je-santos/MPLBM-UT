% PengRobinson.m : calculates the compressibility factor,fugacity coefficient and density 
% of a pure compound with the Peng Robinson equation of state (PR EOS)
% Authors: Piñero.R, Serna.J.G, Martin. A.
%
% function result = PengRobinson(T,P,Tc,Pc,w,MW,Liquido)
% Parameters: T,P,w,Tc,Pc,w,MW,Liquido
% T: Temperature [=] K                                          
% P: Presure [=] Pa                                             
% Tc: critical temperature [=] K                               
% Pc: critical presure [=] Pa                                   
% w: accentic factor
% MW: molar weigth [=] kg/mol
% Liquido:  if Liquido == 1, then calculates liquid fugacity;  
%           if Liquido == 0 then calculates vapor fugacity
% Example:
% [Z fhi density] = PengRobinson(273,2*1.013*1e5,304.21,7.382*1e6,0.225,0.044,1)

%[Z fhi density] = PengRobinson(400,12e6,190.564,4.5992*1e6,0.01142,0.016,1)


function [Z,fhi,density] = PengRobinson(T,P,Tc,Pc,w,MW,Liquido)

R = 8.314; % gas constant [=] J/(mol K)

% Reduced variables
Tr = T/Tc ;
Pr = P/Pc ;

% Parameters of the EOS for a pure component
m = 0.37464 + 1.54226*w - 0.26992*w^2;
alfa = (1 + m*(1 - sqrt(Tr)))^2;
a = 0.45724*(R*Tc)^2/Pc*alfa;
b = 0.0778*R*Tc/Pc;
A = a*P/(R*T)^2;
B = b*P/(R*T);

% Compressibility factor
Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);

ZR = [];
for i = 1:3
   if isreal(Z(i))
   	ZR = [ZR Z(i)];   
   end
end

if Liquido == 1
    Z = min(ZR);   
else
    Z = max(ZR);
end

% Fugacity coefficient
fhi = exp(Z - 1 - log(Z-B) - A/(2*B*sqrt(2))*log((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B)));
if isreal(fhi)
    density=P*MW/(Z*R*T);
    result = [Z fhi density];
else
    'No real solution for "fhi" is available in this phase'
    result=['N/A' 'N/A' 'N/A'];
end
    
% Bibliography: 
%   - ORBEY. H, SANDLER. I; Modeling Vapor-Liquid Equilibria: cubic equations of state and their mixing rules; Cambridge University Press (1998)
%   - Walas,Stanley. M ; Phase Equilibria in Chemical Engineering ; Boston,
%     Butterworth Publishers (1984)                 