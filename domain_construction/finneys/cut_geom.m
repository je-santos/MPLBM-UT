function [new_im, phi] = cut_geom(im, vol_len)
    i_shape   = size(im);
    
    if sum(i_shape==vol_len)==3
        new_im = im;
    else
        first_len = vol_len/2;
        last_len  = vol_len/2-1;

        new_im = im( fix(i_shape(1)/2)-first_len:fix(i_shape(1)/2)+last_len, ...
                     fix(i_shape(2)/2)-first_len:fix(i_shape(2)/2)+last_len, ...
                     fix(i_shape(3)/2)-first_len:fix(i_shape(3)/2)+last_len );
    end
             
    new_im = padarray(new_im,[0,0,1],0,'both'); % add empty voxels for BCs
    [new_im, phi] = eliminate_isolatedRegions(new_im,6); % erase cul-de-sac pores
end
