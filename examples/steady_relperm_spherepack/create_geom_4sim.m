% creates a geometry for simulation

%% opening the file
addpath ('../../pre-processing') %pre-precesing libraries
im_loc = '../../domain_construction/Periodic_spherepacks/';
unzip([im_loc 'periodic_spherepacks_data.zip'], ...
      [im_loc 'periodic_spherepacks_data'])  % unzips the files that contain the geoemtry

d_size = 200; %voxels each side
fp = dlmread([im_loc 'periodic_spherepacks_data/1_200.txt']); %read raw file
fp = reshape(fp, d_size,d_size,d_size);

%% selecting a smaller subset
print_size = 200; %size of the Finneypack subset (in voxels per side)
fp_printing = fp(1:print_size, 1:print_size, 1:print_size);
figure();imagesc(fp_printing(:,:,uint8(print_size/2)));
title('Cross-section of the simulation subset')

%% eliminating non-connected regions 
connect = 6; % pixel connectivity 6, 18, 26
fp = eliminate_isolatedRegions(fp, connect); %for better convergence

%% making a computationally efficent domain for sim
name = ['spheres4Palabos'];
add_mesh   = false;% add a neutral-wet mesh at the end of the domain
num_slices = 0;    % add n empty slices at the beggining and end of domain 
                   % for pressure bcs
swapXZ = true;     % Swap x and z data if needed to ensure Palabos simulation in Z-direction              
scale_2 = false;   % Double the grain (pore) size if needed to prevent single pixel throats
                   % for tight/ low porosity geometries                   

palabos_3Dmat = create_geom_edist(fp_printing,name,num_slices, add_mesh, ...
                                    swapXZ, scale_2);  
                                    %provides a very  efficient 
                                    %geometry for Palabos
                               
%% Mixed Wettability (the user could experiment with this)                                
rng(123)                                    
rnd_array = rand(size(palabos_3Dmat) );

palabos_3Dmat_mixedWet = palabos_3Dmat;
palabos_3Dmat_mixedWet(palabos_3Dmat==1 & rnd_array>0.5)=3;

