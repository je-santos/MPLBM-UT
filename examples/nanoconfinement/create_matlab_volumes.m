% the following script loads a project downloaded from the Digital Rocks
% Portal and prepares it for simulation

% why in matlab? I don't have a good answer for that question

% considerations: I downloaded the projects that interest me and saved them
% in a folder as $project_number_$sample_number
% ie: this sample
% https://www.digitalrocksportal.org/projects/317/origin_data/1354/
% was saved as 317_01.raw
% (.dat extensions were switched to .raw)

% this code might upset some programmers since there are a lot of
% statements that are repeated several times. Each project is slightly
% different and it comes from different people around the world, so I had
% to make sure that everything looked right before simulating the samples
% (which takes many days in many nodes)

% Some things that happen here:
% make sure that the image is labeled as 0 (pore) and 1 (solid)
% crop to 256^3 and 480^3 (when possible)
% and add 2 empty slices in the z-dir
% this is to fit the BCs
% there are some projects involving more than one phase, so I saved each
% phase in a separete file

% last thing,
% here, we will assume that z is connected and that is is the flow dir

global im_save

im_dir  = 'raw_volumes';
im_save = 'binary_volumes_real'; 


for num=339:339 %project number
    
    if num==10
        im_size = 650;
        fb = fopen([im_dir '/' num2str(num) '_01.raw' ],'r');
        im = reshape(fread(fb,im_size^3), im_size,im_size,im_size)/255;
        save_files(im, num, 1)
    end
    
    if num==16
        im_size = 512;
        for i=1:2
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,im_size^3), im_size,im_size,im_size);
            save_files(im, num, i)
        end
    end
    
    
    if num==31
        im_size = [800 220 800];
        for i=1:4
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size)/255;
            im = padarray(im, [0,580/2,0],1,'both');
            save_files(im, num, i)
        end
    end
    
    
    if num==57
        for i=1:7
            im_size = [480 480 480];
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size)/3;
            save_files(im, num, i)
        end
    end
    
    
    
    if num==58
        im_size = [1000 1000 1001];
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size)/9;
        save_files(im, num, 1)
    end
    
    
    if num==65
        im_size = [501 501 501];
        for i=1:5
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size)/255;
            save_files(im, num, i)
        end
    end
    
    if num==69
        im = zeros(255,255,255);
        for i=0:255-1
            im(:,:,i+1) = ~imread( strcat('raw_volumes/69_01/', ...
                                        sprintf('Stack%04d.tif',i ) ) );
        end
        im = padarray(im, [1,1,1],'replicate','pre');
        %im = padarray(im, [1,1,1],'replicate','post'); 
        save_files(im, num, 1)
    end
    
    
    if num==72
        im_size = [1000 1000 1000];
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size);
        im(im==1)=0;
        im(im==9)=1;
        save_files(im, num, 1)
   
        %% subresolution
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size);
        im(im==1)=1;
        im(im==9)=1;
        save_files(im, num, 2)
    end
    
    
    if num==73
        im_size = [501 482 600];
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size);
        im(im==1)=0;
        im(im==9)=1;
        save_files(im, num, 1)
        
        %% subresolution
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size);
        im(im==1)=1;
        im(im==9)=1;
        save_files(im, num, 2)
    end
    
    if num==135
        im_size = [1265 1265 1800];
        for i=0:2
            if i==1
                continue
            end
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size);
            im(im==255)=1;
            im(im==2  )=1;
            save_files(im, num, i)
        end
        
        for i=0:2
            if i==1
                continue
            end
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size);
            im(im==255)=1;
            im(im==0  )=1;
            im(im==2  )=0;
            save_files(im, num, i+4)
        end
        
        i=3;
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size)), im_size);
        im(im==255)=1;
        save_files(im, num, i)
        
    end
    

    if num==172
        im_size = [601 594 1311];
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
        im = ~reshape(fread(fb,prod(im_size),'uint16'), im_size);
        save_files(im, num, 1)
        
        
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(2) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size),'uint16'), im_size);
        im_tmp = im;
        im = (~im) + (im_tmp==2);
        save_files(im, num, 2)
        
        fb = fopen([im_dir '/' num2str(num) '_0' num2str(2) '.raw' ],'r');
        im = reshape(fread(fb,prod(im_size),'uint16'), im_size);
        im_tmp = im;
        im = (~im) + (im_tmp==1);
        save_files(im, num, 3)
    end
    
    if num==204
       packs_dir = 'raw_volumes/204_00/voxelized';
       packs = dir([packs_dir '/*mat']);
       for i=1:numel(packs)
          im = load([packs_dir '/' packs(i).name]);
          im = im.simSpace;
          save_files(im, num, i);
       end
    end
    
    
     if num==207
        im = zeros(861,1024,606);
        for i=1:606
            im(:,:,i+1) = imread('raw_volumes/207_01.tif',i);
        end
        im = ~(im==1)*1.0;
        save_files(im, num, 1)
        
        SE = strel('sphere',1); %mask
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 2)
        
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 3)
        
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 4)
        
    end
    
    
    if num==297
            im_size = [519 557 861];
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size)/255;
            save_files(im, num, 1)
    end

    if num==276
            im_size = [300 300 300];
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(1) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size);
            save_files(im, num, 1)
    end
    
    if num==317
        im_size = [1000 1000 1000];
        for i=1:11
            fb = fopen([im_dir '/' num2str(num) '_0' num2str(i) '.raw' ],'r');
            im = reshape(fread(fb,prod(im_size)), im_size);
            save_files(im, num, i)
        end
    end
    
    if num==339
        for i=1:4
            im = imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'], 1);
            for slice = 2:1094
                im = cat(3, im, ...
                    imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'],slice));
            end
            
            im = im/255;
            save_files(im, num, i)
        end
    end
    
    
    if num==344
        for i=1:4
            im = imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'], 1);
            for slice = 2:1000
                im = cat(3, im, ...
            imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'],slice));
            end
            
            im_tmp = im;
            im = (~im) + (im==1) + (im==5);
            save_files(im, num, i)
        end
        
        
        
        for i=1:4
            im = imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'], 1);
            for slice = 2:1000
                im = cat(3, im, ...
            imread( [im_dir '/' num2str(num) '_0' num2str(i) '.tif'],slice));
            end
            im_tmp = im;
            im = (~im) + (im==2) + (im==5);
            save_files(im, num, i+4)
        end
        
    end
    
    if num==370
        
        zsize = 848; % 850, 848
        im = zeros(919,990,zsize);
        for i=1:zsize
            im(:,:,i+1) = imread('raw_volumes/370_04.tif',i);
        end
        im = (im==0)*1.0;
        save_files(im, num, 1)
        SE = strel('sphere',1); %mask
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 1)
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 1)
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 1)
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 1)
        im = uint8(1-imdilate(~im,SE)); %pore space dilation
        save_files(im, num, 1)
        
        
        
        
        im = zeros(855,934,848);
        for i=1:848
            im(:,:,i+1) = imread('raw_volumes/370_02.tif',i);
        end
        im = (im==0)*1.0;
        save_files(im, num, 2)
        
        
        
        
        
%         
%         
%         SE = strel('sphere',1); %mask
%         im = uint8(1-imdilate(~im,SE)); %pore space dilation
%         save_files(im, num, 2)
%         
%         im = uint8(1-imdilate(~im,SE)); %pore space dilation
%         save_files(im, num, 3)
%         
%         im = uint8(1-imdilate(~im,SE)); %pore space dilation
%         save_files(im, num, 4)
        
    end
    
    
    
end









