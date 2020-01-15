%% Input for the function
dir_name = '1e-2';
domain_size = [75,72,215];

mesh_added = true; % was add a neutral-wet mesh at the end of the domain?
num_slices = 2;    % how many n empty slices at the beggining and end of domain 

%% function
breakthrough = 0;

addpath ('../../pre-processing') %pre-precesing libraries
all_rhos = dir( [dir_name '/rho*dat'] );
all_vels = dir( [dir_name '/vel*dat'] );

domain_loc = dir('input/*dat');
domain     = load( ['input/' domain_loc.name ] );
domain     = reshape(domain, flip(domain_size));

domain = erase_voxels(domain,mesh_added,num_slices);

vol_total = sum(~domain(:)); %total volume available for flow


for pc_inc= 1:size(all_rhos,1)
    
    
   rho = load( [dir_name '/' all_rhos(pc_inc).name] );
   rho = reshape(rho, flip(domain_size));
   
   rho = erase_voxels(rho,mesh_added,num_slices);
   
   fluid_mask   = rho>=1;
   saturation   = sum(fluid_mask(:))/vol_total
   
   if breakthrough == 0
    breakthrough = sum(sum( rho(:,:,end) )) > 1
   end
   
   %figure;imagesc(squeeze(rho(:,int8(domain_size(2)),:)))
   
   vel = load( [dir_name '/' all_vels(pc_inc).name] );
   vel = reshape( vel, [3 flip(domain_size)] );
   vz = squeeze( vel(1,:,:,:) ); %grab the first component of velocity (z)
   clear vel
   
   vz = erase_voxels(vz,mesh_added,num_slices);
   
   mean_vz  = mean( vz(fluid_mask)  )
   mean_rho = mean( rho(fluid_mask) )
   ca       = 0.166667*mean_rho*mean_vz/( 0.15*0.8 )
   
   % final domain for relative permeability
   domain4kr_f1 = ~fluid_mask;
   domain4kr_f2 = ~domain + ~fluid_mask;
   
   name='tmp';
   
   %palabos_kr_f1 = create_geom_edist(domain4kr_f1,name,0,false,false,false);  
   
                                
                    %print a file with headers with runnum, Pc, S, Bk, Ca, Tortuosity
   
end



function array = erase_voxels(array,mesh_added,num_slices)
   array(:,:,1:num_slices                       )=[];  %removes the first slices
   array(:,:,end-num_slices-(mesh_added*2)+1:end)=[];  %removes the last slices + mesh +1 more
end

