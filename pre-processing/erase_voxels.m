function array = erase_voxels(array,mesh_added,num_slices)
    array(:,:,1:num_slices                       )=[];  %removes the first slices
    array(:,:,end-num_slices-(mesh_added*1)+1:end)=[];  %removes the last slices + mesh +1 more
end