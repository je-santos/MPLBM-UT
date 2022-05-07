function []=crop_save(crop_size, im, num, i)
    i_shape = size(im);
    if sum( i_shape >= crop_size) == 3
            [new_im, phi] = cut_geom(im, crop_size);
            if phi>0.01
                write_im(new_im, num, i)
            else
                disp('The porosity is lower than 1%')
            end
    else
        disp(['The geometry is too small to crop ' num2str(crop_size)])
    end
end