%% Experiment name and options
create_batch=true; %if true creates batchjob

send_batch=false; %if true executes batchjob, otherwise it runs it in the node


exp_name='imb_1e-5';



%% Input Data
%dp=0.13; %pressure change (Inlet-Outlet) prev 0.4/0.1 ->0.09/0.07
dp=0.0
output_folder=exp_name;
system(['mkdir ',output_folder]);
fileID = fopen([output_folder '/input.txt'],'w');
path_output=[output_folder,'/ '];
omega1= '1.0 ';
omega2= '1.0 ';
nx='106 '; %depth
ny='92 ';
nz='254 ';
g='0.9 ';
force_f1='0.0 ';
force_f2='0.0 ';
g_ads_f1_s1='0.4 '; % the f2 is given by its negative, (0.4nw,-0.4w)
g_ads_f1_s2='0.0 ';
domain='254x92x106.dat ';
max_it='1000000 ';
save_it='5000 ';
vtk_it='1000000 ';
info_it='500 ';%displays info and checks for convergence
convergence='0.00001 ';
start_no='0 '; %no of first realization
run_no='200 '; %number of densities
nx1f1='1 '; %initial fluid1 distribution
nx2f1='2 '; 
ny1f1='0 '; 
ny2f1='92 '; 
nz1f1='0 '; 
nz2f1='254 '; 
nx1f2='3 '; %initial fluid2 distribution
nx2f2='105 '; 
ny1f2='0 '; 
ny2f2='92 '; 
nz1f2='0 '; 
nz2f2='254 ';
rho_no_fluid='0.06 '; %check if works for static case
cap_factor='1.0 ';

%Writting Pressure
runs=str2num(run_no);
pressure=[];
pressure(1)=1.5;
pres_str(1)={num2str(pressure(1))};
pres_str2(1)={num2str(2.00)};
for i=1:(runs-1)
   pressure(i+1)=pressure(i)+(dp-(pressure(1)-2))/(runs-1);
   pres_str(i+1)={num2str(pressure(i+1))};
   pres_str2(i+1)={num2str(2.0,'%f')};
end
pres_str1='1.5'; %initial density
pres_str22='2.0';
for i=1:runs-1
   pres_str1=strcat(pres_str1,{' '},pres_str(i+1));
   pres_str22=strcat(pres_str22,{' '},pres_str2(i+1));
end
pres_str1=strcat(pres_str1,{' '});
%pres_str22=strcat(pres_str22,{' '});
rho1=cell2mat(pres_str1);
rho2=cell2mat(pres_str22);   %Uncomment for more than 1 density

%rho1='2.04 ';
%rho2='2.0';

q=char(39); %single quote for printing

%% Executing string
%str_sc=['system([',q,'ibrun ./ShanChen ',q,' ',path_output,omega1,omega2,nx,ny,nz,g,...
%    force_f1,force_f2,g_ads_f1_s1,g_ads_f1_s2,domain,max_it,save_it,vtk_it,info_it,...
%    convergence,start_no,run_no,nx1f1,nx2f1,ny1f1,ny2f1,nz1f1,nz2f1,nx1f2,nx2f2,ny1f2,ny2f2,nz1f2,nz2f2,...
%    rho_no_fluid,cap_factor,rho1,rho2, '])']

str_sc=['ibrun ./ShanChen ',path_output,omega1,omega2,nx,ny,nz,g,...
    force_f1,force_f2,g_ads_f1_s1,g_ads_f1_s2,domain,max_it,save_it,vtk_it,info_it,...
    convergence,start_no,run_no,nx1f1,nx2f1,ny1f1,ny2f1,nz1f1,nz2f1,nx1f2,nx2f2,ny1f2,ny2f2,nz1f2,nz2f2,...
    rho_no_fluid,cap_factor,rho1,rho2];

fprintf(fileID,str_sc);

%% Batch
if create_batch==true
    fileID_batch = fopen(['batchjob_',exp_name],'w');
    fprintf(fileID_batch,'#!/bin/bash \n#SBATCH -J myjob');
    fprintf(fileID_batch,'          # job name \n');
    fprintf(fileID_batch,'#SBATCH -o ');
    fprintf(fileID_batch,output_folder);
    fprintf(fileID_batch,'%s','/myMPI.o%j        # output and error file name (%j expands to jobID)');
    fprintf(fileID_batch,         '\n');
    fprintf(fileID_batch,'#SBATCH -N 10                # number of nodes requested \n');
    fprintf(fileID_batch,'#SBATCH -n 680               # total number of mpi tasks requested \n');
    fprintf(fileID_batch,'#SBATCH -p normal      # queue (partition) -- normal, development, etc. \n');
    fprintf(fileID_batch,'#SBATCH -t 24:00:00         # run time (hh:mm:ss) - 1.5 hours \n\n');
    fprintf(fileID_batch,'# Slurm email notifications are now working on Lonestar 5 \n');
    fprintf(fileID_batch,'#SBATCH --mail-user=jesantos@utexas.edu \n');
    fprintf(fileID_batch,'#SBATCH --mail-type=end     # email me when the job finishes \n');
    fprintf(fileID_batch,'# run the executable named a.out \n');
    fprintf(fileID_batch,str_sc);
    fclose(fileID_batch);
    disp('Batchjob created')
end
fclose(fileID);

if send_batch==true
    system(['sbatch -A pge- batchjob_',exp_name])
else
    system('idev') %just in case
    system(str_sc)
end

