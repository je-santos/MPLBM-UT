function [] = writeSlurmJob(geom_name, geom_loc, pressure, sim_size, saveto)

%======================================================================
%input parameters
%======================================================================
%clear
cores_per_node = 36;

p0_P=pressure*1.0e6;        %pressure, [Pa]
T0_P=400.0;                 %temperature, [k]
dx=0.5e-9;                  %lattice resolution, [m/lu]

H_L=sim_size-1;             %characterisitic length in lattice, [lu]   %changefordiffsize:)
H_P=dx*H_L;                 %characteristic length in physical unit, [m]
%=======================================================================
pressgra_P=0.1E6;           %pressure gradient, [Pa/m]
%-----------------------------------------------------------------------

%%------------ physical unit and lattice unit --------------------------
Rig=8.314;                  %ideal gas constant, [J/(mol·K)]
d_CH4=0.375e-9;             %diameter of methane molecule
M_CH4=16.04e-3;             %molar weigth, [kg/mol]
NA=6.02e23;                 %Avogadro’s number,[1/mol]
mass_CH4=M_CH4/NA;          %mass of methane molecule
w=0.01142;                  %accentic factor

%------------ critical parameters for methane --------------------------
Tcr_P=190.564;              %critical temperature in physical unit,NIST, [K]
rhocr_P=162.66;             %critical density in physical unit,NIST,     [kg/m^3]
pcr_P=4.5992E6;             %critical pressure in physical unit,NIST,    [Pa]

%=======================================================================

[Z,fhi,rho0_P] = PengRobinson(T0_P,p0_P,Tcr_P,pcr_P,w,M_CH4,0); %a function to obtain density, based on P-R EOS
Num=rho0_P/mass_CH4;                %molecular number density of methane,

%dynamic viscousity from Lee et al.1966, [Pa*s]
K_Lee=(9.379+0.01607*M_CH4)*(T0_P^1.5)/(209.2+19.26*M_CH4+T0_P);
X_Lee=3.448+986.4/T0_P+0.01009*M_CH4;
Y_Lee=2.447-0.2224*X_Lee;
mu_P=1.0d-7*K_Lee*exp(X_Lee*((rho0_P/1000.0)^Y_Lee));

%kinetic viscosity
nu_P=mu_P/rho0_P;                    %Kinematic viscosity, [m^2/s] 

rho0_L=1.0;                          %free gas density in lattice, [lu]  
rho0_PL=rho0_P/rho0_L;               %density scale

H_PL=H_P/H_L;                        %length scale, [m/lu]
MFP_P=1.0/(sqrt(2.0)*Num*pi*d_CH4^2);%mean free path in physical unit, [m]

MFP_L=MFP_P/H_PL;                    %mean free path in lattice unit, [lu]	
Kn=MFP_P/H_P;                        %Knudsen number

nu_L=sqrt(2.0/pi/3.0)*MFP_L;         %Kinetic viscosity in lattice unit,[lu] 
nu_PL=nu_P/nu_L;                     %Kinetic viscosity scale, [(m^2/m)/lu]

usurf_PL=nu_PL/H_PL;                 %velocity scale, [(m/s)/lu]
usurf_P=0.0;
usurf_L=usurf_P/usurf_PL;            %velocity in lattice unit, [lu]

accelz_P=pressgra_P/rho0_P;          %accelation in z direction, [m/s^2]
accelz_PL=usurf_PL^2/H_PL;           %accelation scale, [(m/s^2)/lu]
accelz_L=accelz_P/accelz_PL;         %accelation in lattice unit, [lu]

%%------------------------------------------------------------------------

jobName=[geom_name '_' num2str(pressure)];
batchName=['JOB_' jobName];
numNode=4*2;
tHour=16; 
%taccAcct='pge-fracture'; 
mpiExec='lbm_lev';

nameSim=[ saveto '/' jobName ];  %This is the name of the output folder
mkdir(nameSim)

typeConv='vel'; %currently two choices either 'vel' or 'rho' readCVolout_Reconstruct
NumPx='4';
NumPy='6';
%NumPz='6'; 
NumPz='12'; 

totalprocs = num2str( str2num(NumPx)*str2num(NumPy)*str2num(NumPz) );

geomFile= geom_loc; %Ensure this matches lattice file
lx=num2str(sim_size); 
ly=num2str(sim_size); 
lz=sim_size+2; 
perX_='0'; 
perY_='0'; 
perZ_='1';
accelx='0.0';
accely='0.0';
accelz='0.00000001';%num2str(accelz_L); 

gradP=0.0;
Px_='0'; 
Px_LinP_='0';
PxIn='0.0'; 
PxOut='0.0';
Py_='0'; 
Py_LinP_='0';
PyIn='0.0';
PyOut='0.0'; 

Pz_='0'; 
Pz_LinP_='0';
PzIn=1/3+(gradP*lz); %density = pressure*3
PzOut=1/3;

resBx_='0'; 
resRhoxIn='0.0'; 
resRhoxOut='0.0';
resBy_='0'; 
resRhoyIn='0.0';
resRhoyOut='0.0'; 
resBz_='0'; 
resRhozIn='0.0';    %match to input rho file
resRhozOut='0.0';
isoNum='1';         %1 for 4th order isotropy, 2 for 8th order isotropy, 3 for 10th order isotropy
globalVisc_='0';    %set globalVisc_ to 1 to turn off local effective viscosity
nu='0.1666666666';
mfp=num2str(MFP_L); %in lattice units, used to calculate the unbound viscosity
stdBB_='0'; %Turn on standard BB, this is unlikely to be done for nanoconfinement simulations, results in no-slip boundary
uzWall=num2str(usurf_L); %surface diffusion velocity in lattice unit
%uzWall='0.0'if ignoring the adsorbed gas  
rBB='0.8909';   %Following Chai et al. 2010,  rSuga(sliceCount)=rBB
%%---------------------------------------------
GSC='0'; %%Gf=Gads, GSC=0 in Shan-Chen if ignoring the adsorbed gas. But GSC~=0 in P-R, even if GSC doesn't work in P-R EOS
rhoFluid=num2str(rho0_L);   %fluid density in lattice unit
rhoFluidran='0.000';  
rhoWall=num2str(rho0_L); %rhoWall=rhoFluid if ignoring the adsorbed gas

%------------------------------------------------------
readRhoIn_='0'; %1 on, 0 off
rhoInName='dummy';
convCrit='0.00001'; 
maxT='1000001'; %max number iterations
convT='1000'; %Check for convergence
outT='500000'; %outputs at mod(outT) 
ux_='1'; %1 to output the ux velocity and so on...
uy_='1'; 
uz_='1'; 
rho_='1'; 
fIn_='0'; 
fOut_='0'; 
geom_='0';
%%%%%%%%%%%%%%%%%%%%%%%%%%%End User Inputs%%%%%%%%%%%%%%%%

% mkdir([nameSim]);
fid = fopen([ saveto '/' batchName], 'w');    % open the output file to write in
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '\n');
fprintf(fid, ['#SBATCH -J ' jobName '\n']);
fprintf(fid, ['#SBATCH -o ' jobName '.o%%j\n']);
fprintf(fid, ['#SBATCH -e ' jobName '.er\n']);
fprintf(fid, '#SBATCH -N %i\n', numNode);
fprintf(fid, '#SBATCH -n %i\n', (numNode*cores_per_node));
%fprintf(fid, '#SBATCH -p normal\n');
hh=num2str(floor(tHour),'%0.2i');
mm=num2str((tHour-floor(tHour))*60, '%0.2i');
ss='00';
tim=[hh ':' mm ':' ss];
fprintf(fid, ['#SBATCH -t ' tim '\n']);
fprintf(fid, '\n');

fprintf(fid, 'module purge \n');
fprintf(fid, 'module load cmake/3.16.2 \n');
fprintf(fid, 'module load intel-mpi/2019.9.304 \n');
fprintf(fid, 'module load intel/15.0.5 \n');

fprintf(fid, '\n');
%fprintf(fid, '\n');
%fprintf(fid, ' set -x\n');
%fprintf(fid, '\n');

tmp = split( nameSim, '/');
nameSim = tmp{end};

fprintf(fid, ['mpirun -n ' totalprocs ' ./' mpiExec ' ' nameSim ' ' geomFile ' ' NumPx ' ' NumPy ' ' NumPz ' ' lx ' ' ly ' ' num2str(lz) ' ' perX_ ' ' perY_ ' ' perZ_ ' ' ...
                accelx ' ' accely ' ' accelz ' ' Px_ ' ' Px_LinP_ ' ' PxIn ' ' PxOut ' ' Py_ ' ' Py_LinP_ ' ' PyIn ' ' PyOut ' ' Pz_ ' ' Pz_LinP_ ' ' num2str(PzIn,'%.10f') ' ' num2str(PzOut,'%.10f') ' ' ...
                resBx_ ' ' resRhoxIn ' ' resRhoxOut ' ' resBy_ ' ' resRhoyIn ' ' resRhoyOut ' ' resBz_ ' ' resRhozIn ' ' resRhozOut ' ' ...
                isoNum ' ' globalVisc_ ' ' nu ' ' mfp ' ' stdBB_ ' ' uzWall ' ' rBB ' ' GSC ' ' ...
                rhoFluid ' ' rhoFluidran ' ' rhoWall ' ' readRhoIn_ ' ' rhoInName ' ' convCrit ' ' typeConv ' ' maxT ' ' convT ' ' outT ' ' ...
                ux_ ' ' uy_ ' ' uz_ ' ' rho_ ' ' fIn_ ' ' fOut_ ' ' geom_]);

fclose(fid);


end
