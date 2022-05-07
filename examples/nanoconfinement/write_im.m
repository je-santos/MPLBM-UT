function [] = write_im(new_im, num, i)
    global im_save
    i_shape = size(new_im);
    x_size  = i_shape(1);
    
    bin = uint8(new_im);
    
   
    save([im_save '/' num2str(num) '_0' num2str(i) '_' num2str(x_size)], ...
             'bin')
    
    figure;imagesc(new_im(:,:,100));pause(.1);
    
end