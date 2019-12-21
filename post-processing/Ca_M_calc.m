%Calculating viscosity ratios (M) and Capillary numbers (Ca) from
%simulation results


%User inputs:

absK= 0.1; % Value from 1 phase simulation
sigma= 0.15; % Interfacial tension
%Update path to folder where output .dat files from 2-phase simulation are stored
directory='C:\Users\Abhishek\Desktop\shanchen_updated_try\MultiphasePorousMediaPalabos\examples\relative_permeability_spherepack\tmp\'; 

%%%%%%%%%%%%%%%
fileID = fopen([directory 'output.dat']);
C = textscan(fileID,'%s','Delimiter','=');
C1=C{1};
C2 = str2double(C1);

Index_p = find(contains(C1,'Pressure difference'));
pressure_num=0;

for i=1:numel(Index_p)
    pressure_num=pressure_num+1;
    current_index=Index_p(i);
    pressure=C2(current_index+1);
    delP(i,1)=pressure;
end

Index_d = find(contains(C1,'Inlet density'));
rho1=C2(Index_d+1);  %inlet density

Index_d1 = find(contains(C1,'Dissolved density'));
dis_rho=C2(Index_d1+1);

Index_gads = find(contains(C1,'Gads_f1_s1'));
Gads_f1_s1=C2(Index_gads+1);
Gads_f1_s2=-Gads_f1_s1;

Index_gc = find(contains(C1,'Gc'));
Gc=C2(Index_gc+1);

Index_l = find(contains(C1,'Geometry flow length'));
flow_length=C2(Index_l+1);

Index_v1 = find(contains(C1,'Kinematic viscosity f1'));
kine_visco1=C2(Index_v1+1);

Index_v2 = find(contains(C1,'Kinematic viscosity f2'));
kine_visco2=C2(Index_v2+1);

delRho = delP*3; % Density difference
rho2=rho1-delRho;% Density of fluid 2 (outlet)

dyna_vis_1= kine_visco1*rho1;
dyna_vis_2= kine_visco2*rho2;

%Calculating velocity of fluid 1 by Darcy's law
vel_1 = absK*delP/(dyna_vis_1*flow_length);

%Calculating contact angle according to Huang et al. (2007)
cosTheta = abs((Gads_f1_s2-Gads_f1_s1)/(Gc*(rho1-dis_rho)*0.5));


M=dyna_vis_1./dyna_vis_2; %Viscosity ratio
Ca=(abs(dyna_vis_1*vel_1)/(sigma*cosTheta));  %Capillary number


CaM=cat(1,M,Ca);
csvwrite([directory '\' 'Ca_M' '.csv'],CaM); % Saves values for each pressure increment
