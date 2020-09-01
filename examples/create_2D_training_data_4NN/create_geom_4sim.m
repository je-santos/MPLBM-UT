% creates geometries for simulation

%% opening the images
addpath ('../../pre-processing') %pre-precesing libraries
geom_loc = '../../domain_construction/micromodels/';
mkdir('input');
% opening the micromodels
berea   = Tiff([geom_loc 'Berea.tif']); berea = logical(read(berea))*1;
cement  = imread([geom_loc 'cementation_mask.PNG']);cement=logical(cement(:,:,1))*1;
disol   = imread([geom_loc 'dissolution_mask.PNG']);disol=logical(disol(:,:,1))*1;
rock    = imread([geom_loc 'Fracture_mask.PNG']);rock=logical(rock(:,:,1))*1;
frac    = Tiff([geom_loc 'Fracture_mask2.tif']);frac=~read(frac)*1;
snake   = Tiff([geom_loc 'Snake_EtchedPorosity.tif']);snake=read(snake);snake=snake(:,:,1);
snake2  = imread([geom_loc 'SnakePorosity.PNG']);snake2=(snake2(:,:,3)>100)*1;
vug     = imread([geom_loc 'vug_mask.PNG']);vug=logical(vug(:,:,1))*1;

%% arrays of images for subsecuent operations
for_downsizing = {cement, disol, frac}; 
for_resampling = {berea, cement, disol, rock, frac, snake, snake2, vug};
for_resizing   = {berea, snake, snake2};
for_sim        = {};

%% further image processing
rng(123123)
print_size_x = 128; 
print_size_y = 256;

%% rescaling the images
for im=1:size(for_downsizing,2)
    for_resampling{end+1} = round(imresize(for_downsizing{im}, 1/2));
    for_resampling{end+1} = round(imresize(for_downsizing{im}, 1/4));
end

for im=1:size(for_resizing,2)
    for_resampling{end+1} = round(imresize(for_resizing{im}, 1.5));
    %for_resampling{end+1} = round(imresize(for_resizing{im}, 4));
end

for im=1:size(for_resizing,2)
    for_sim{end+1} = imresize(for_resizing{im}, [print_size_x,print_size_y],'nearest');
    for_sim{end+1} = imresize(for_resizing{im}, [print_size_x,print_size_y],'nearest');
end


%% eroding and dilating the images
SE = strel('sphere',1); %mask
for im=1:size(for_resampling,2)
    tmp_im  =  for_resampling{im};
    tmp_im2 = ~for_resampling{im};
    
    for i=1:3
        tmp_im = imdilate(tmp_im,SE); %pore space erosion
    end
    for_resampling{end+1} = tmp_im;
    
   
    for i=1:3
        tmp_im = imdilate(tmp_im,SE); %pore space erosion
    end
    for_resampling{end+1} = tmp_im;
    
    for i=1:3
        tmp_im = imdilate(tmp_im,SE); %pore space erosion
    end
    for_resampling{end+1} = tmp_im;
    
    % dilating the grains
    for i=1:2
        tmp_im2 = imdilate(tmp_im2,SE); %pore space erosion
    end
    for_resampling{end+1} = ~tmp_im2;
    
    for i=1:2
        tmp_im2 = imdilate(tmp_im2,SE); %pore space erosion
    end
    for_resampling{end+1} = ~tmp_im2;
 
end


%% create domains for simulation
num_resamples = 10;

for im=1:size(for_resampling,2)
    im_size =  size(for_resampling{im});
    for resamples=1:num_resamples
        
        if im_size(1)-print_size_x>0 && im_size(2)-print_size_y>0
            xi = randi([1, im_size(1)-print_size_x],1);
            yi = randi([1, im_size(2)-print_size_y],1);
            tmp_im = for_resampling{im};
            for_sim{end+1} = tmp_im(xi:xi+print_size_x-1,yi:yi+print_size_y-1);
        end
        
        if im_size(1)-print_size_y>0 && im_size(2)-print_size_x>0
            xi = randi([1, im_size(2)-print_size_x],1);
            yi = randi([1, im_size(1)-print_size_y],1);
            tmp_im = for_resampling{im}';
            for_sim{end+1} = tmp_im(xi:xi+print_size_x-1,yi:yi+print_size_y-1);
        end
        
        
    end
    
end


% create LBM domains
geom.name = ['data_4NN'];
geom.print_size = true;
geom.add_mesh   = false; % add a neutral-wet mesh at the end of the domain
geom.num_slices = 0;     % add n empty slices at the beggining and end of domain 
                         % for pressure bcs
geom.swapXZ = false;     % Swap x and z data if needed to ensure Palabos simulation in Z-direction              
geom.scale_2 = false;    % Double the grain (pore) size if needed to prevent single pixel throats
                         % for tight/ low porosity geometries    

sim_domains={};
for im=1:size(for_sim,2)
    phi = sum(for_sim{im}>0,'all')/numel(for_sim{im})*100;
    if phi > 10 && phi<75
        tmp_im = repmat(~for_sim{im},1,1,3);
        
        blank_slice  = tmp_im(:,1:1,:)*0 ;
        tmp_im = cat(2, blank_slice, tmp_im);
        tmp_im = cat(2, tmp_im, blank_slice);
        
        sim_domains{end+1} = create_geom_edist(tmp_im, geom); %provides a very  efficient geometry 4sim
        
        tmp = sim_domains{end};
        imagesc(tmp(:,:,2))
        title(['Im ' num2str(im)])
        pause(.5)
    end
end

