%% Input for the function
dir_name = '1e-2';
domain_size = [75,72,215];

mesh_added = true; % was add a neutral-wet mesh at the end of the domain?
num_slices = 2;    % how many n empty slices at the beggining and end of domain 

%% function
addpath ('../../pre-processing') %pre-precesing libraries
all_rhos = dir( [dir_name '/rho*dat'] );
all_vels = dir( [dir_name '/vel*dat'] );


for pc_inc= 1:size(all_rhos,1)*0+3
    
   domain_loc = dir('input/*dat');
   domain     = load( ['input/' domain_loc.name ] );
   domain     = reshape(domain, flip(domain_size));
   
   domain = erase_voxels(domain,mesh_added,num_slices);
    
   rho = load( [dir_name '/' all_rhos(pc_inc).name] );
   rho = reshape(rho, flip(domain_size));
   
   rho = erase_voxels(rho,mesh_added,num_slices);
   
   fluid_mask = rho>1;
   
   vel = load( [dir_name '/' all_vels(pc_inc).name] );
   vel = reshape( vel, [3 flip(domain_size)] );
   vz = squeeze( vel(1,:,:,:) ); %grab the first component of velocity (z)
   clear vel
   
   vz = erase_voxels(vz,mesh_added,num_slices);
   
   mean_vz = mean(vz(fluid_mask));
   
   % final domain for relative permeability
   domain4kr = fluid_mask;
    
end




function array = erase_voxels(array,mesh_added,num_slices)
   array(:,:,1:num_slices                       )=[];  %removes the first slices
   array(:,:,end-num_slices-(mesh_added*2)+1:end)=[];  %removes the last slices + mesh +1 more
end

