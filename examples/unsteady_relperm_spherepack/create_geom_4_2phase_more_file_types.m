% This creates a geometry for simulation from different types-
% Type 1) Raw image file
% Type 2) MATLAB (.MAT) file
% Type 3) image slices
% Type 4) sphere centre coordinates and radii

addpath ('../../pre-processing') %pre-procesing libraries

%% Define file type (1/2/3/4)
file_type = 1;
d_size = 500; % voxels each side (specify for raw file)

%% opening the file

% raw image geometry file
if file_type == 1   
    f1 = fopen('input/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw','r'); %read raw file
    fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
    fp = reshape(fp, d_size,d_size,d_size);

% .MAT geometry file    
elseif file_type == 2 
    fp=load('input/*.mat'); %replace * with file name if known

% image slices file by type    
elseif file_type == 3
    fp = read_img_slices('input/*.tiff'); %replace image slice file type eg(.png | .tiff | .jpg)
    
% sphere centre coordinates and radii    
elseif file_type == 4
    res = 50;  % Resolution for scaling coordinates to voxels
    fp = read_csv_file('input/*.csv', res); %replace * with file name if known
    
else
    print('Invalid file_type. Please choose from 1-4');
end
    
%% selecting a smaller subset
print_size = 200; %size of the Finneypack subset (in voxels per side)
fp_printing = fp(1:print_size, 1:print_size, 1:print_size);
figure();imagesc(fp_printing(:,:,uint8(print_size/2)));
title('Cross-section of the simulation subset')

%% eliminating non-connected regions 
connect = 6; % pixel connectivity 6, 18, 26
fp_printing = eliminate_isolatedRegions(fp_printing, connect); %for better convergence

%% making a computationally efficent domain for sim
geom.name       = ['spheres4Palabos'];
geom.print_size = true;
geom.add_mesh   = true; % add a neutral-wet mesh at the end of the domain
geom.num_slices = 4;    % add n empty slices at the beggining and end of domain 
                   % for pressure bcs
geom.swapXZ     = true;     % Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2    = false;   % Double the grain (pore) size if needed to prevent single pixel throats
                   % for tight/ low porosity geometries                   

palabos_3Dmat   = create_geom_edist(fp_printing,geom);  
                                    %provides a very  efficient 
                                    %geometry for Palabos


% old version
%palabos_3Dmat = mat2dat_4lbm(fp_printing,name,1); %although this function is slow, it 
                                    

%% Mixed Wettability (the user could experiment with this)                                
rng(123)                                    
rnd_array = rand(size(palabos_3Dmat) );

palabos_3Dmat_mixedWet = palabos_3Dmat;
palabos_3Dmat_mixedWet(palabos_3Dmat==1 & rnd_array>0.5)=3;

