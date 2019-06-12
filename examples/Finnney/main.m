addpath ('../../pre-processing') %pre-precesing libraries
d_size = 500;
f1 = fopen('spheres_a10_dx0.04_n500_segmented_unsigned_char.raw','r'); %read raw file
fp = fread(f1, d_size*d_size*d_size,'uint8=>uint8');
fp = reshape(fp, d_size,d_size,d_size);
fp = eliminate_isolatedRegions(fp, connect);
fp = fp(0:100,0:100, 0:100); %subsamples the geometry into 100^3

%% print for palabos

% call a function to get the geom printed in this folder
% i.e.
% printPalabos(fp,'nameForTheFile')
% if you can get an imageSC of the process, that would be great :)