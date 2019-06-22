%% creates a geometry for simulation

addpath ('../../pre-processing') %pre-precesing libraries
d_size = 500; %voxels each side
f1 = fopen('input/spheres_a10_dx0.04_n500_segmented_unsigned_char.raw','r'); %read raw file
fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
fp = reshape(fp, d_size,d_size,d_size);

connect = 6; % pixel connectivity 6, 18, 26
fp = eliminate_isolatedRegions(fp, connect); %for better convergence

%% print for palabos
name='spheres4Palabos';
[status] = mat2dat_4lbm(fp,name,1); %although this function is slow, it 
                                    %provides a very computationally efficient 
                                    %geometry for Palabos
                                    
% TODO: update faster but less efficient printing function
