% creates a geometry for simulation

%% opening the file
addpath ('../../pre-processing') %pre-precesing libraries
d_size = 480; %voxels each side
f1 = fopen('input/grid_25.bin','r'); %read raw file
fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
fp = reshape(fp, d_size,d_size,d_size);

%% selecting a smaller subset
print_size = 250; %size of the subset (in voxels per side)
fp_printing = fp(1:print_size, 1:print_size, 1:print_size)/max(fp(:));
figure();imagesc(fp_printing(:,:,uint8(print_size/2)));
title('Cross-section of the simulation subset')

%% eliminating non-connected regions 
connect = 6; % pixel connectivity 6, 18, 26
fp_printing = eliminate_isolatedRegions(fp_printing, connect); %for better convergence

%% making a computationally efficent domain for sim
geom.name       = ['carbonate4Palabos'];
geom.print_size = true;
geom.add_mesh   = false;% add a neutral-wet mesh at the end of the domain
geom.num_slices = 1;    % add n empty slices at the beggining and end of domain 
                   % for pressure bcs
geom.swapXZ     = true;     % Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2    = false;   % Double the grain (pore) size if needed to prevent single pixel throats
                   % for tight/ low porosity geometries                   

palabos_3Dmat = create_geom_edist(fp_printing,geom);  
                                    %provides a very  efficient 
                                    %geometry for Palabos


