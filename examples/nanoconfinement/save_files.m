function [] = save_files(im, num, i)
    disp(['Sample ' num2str(num) ' file ' num2str(i)])
    
    if numel(unique(im))>2
        error('The image is not binary')
    end
    
   
    crop_size = 256;
    crop_save(crop_size, im, num, i);
    crop_size = 480;
    crop_save(crop_size, im, num, i);  
end