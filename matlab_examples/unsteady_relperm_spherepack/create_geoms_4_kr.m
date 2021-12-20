addpath ('../../post-processing') %post-processing libraries

%% Input for the function
kr.domain_size  = [200,200,200];
kr.mesh_added   = true; % was a neutral-wet mesh added at the end of the domain?
kr.num_slices   = 4;    % how many n empty slices at the begining and end of domain
kr.input_dir    = 'input';
kr.input_geom   = 'spheres4Palabos';  %This is original geometry in input folder. Output geometries from 2-phase simulation should be put in another folder (eg. tmp)
kr.output_file  = 'info_4kr.txt';
kr.pressure_bcs = true;
% if true, the following info is needed for the Pc-S curve
kr.rho_init_f1  = 0; %initial density of fluid 1 at the boundary
kr.rho_final_f1 = 0; %final density of fluid 1 at the boundary
kr.rho_steps    = 0; %number of pressure increments

domains_4_kr
