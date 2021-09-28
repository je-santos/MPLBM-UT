global im_save

im_dir  = 'synthetic_volumes';
im_save = 'binary_volumes_syn';

%% Finneypack operations
num=1;
f1=fopen('spheres_a10_dx0.04_n500_segmented_unsigned_char.raw','r');
d_size = 500;
im=fread(f1,d_size^3,'uint8=>uint8');
im=reshape(im,d_size,d_size,d_size);


SE = strel('sphere',1); %mask
im_dilated = uint8(1-imdilate(~im,SE)); %pore space dilation
save_files(im_dilated, num, 0)

for i=1:8
    save_files(im, num, i)
    im = imdilate(im,SE); %pore space erosion
end
