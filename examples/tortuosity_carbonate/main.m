%% creates a geometry for simulation

addpath ('../../pre-processing') %pre-precesing libraries
d_size = 480; %voxels each side
f1 = fopen('input/grid_25.bin','r'); %read raw file
fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
fp = reshape(fp, d_size,d_size,d_size);

connect = 6; % pixel connectivity 6, 18, 26
fp = eliminate_isolatedRegions(fp, connect); %for better convergence

%% print for palabos

print_size = 480; %size of the Finneypack subset (in voxels per side)
fp_printing = fp(1:print_size, 1:print_size, 1:print_size);

figure();imagesc(fp_printing(:,:,uint8(print_size/2)));
title('Cross-section of the simulation subset')

name = ['carbonate4Palabos'];

palabos_3Dmat = mat2dat_4lbm(fp_printing,name,1); %although this function is slow, it 
                                    %provides a very computationally efficient 
                                    %geometry for Palabos

