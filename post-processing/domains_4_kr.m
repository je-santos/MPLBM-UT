%% function (kr.)
addpath ('../../pre-processing') %pre-precesing libraries
fileID = fopen(kr.output_file,'w');
breakthrough = 0; %initialize
% look for the subfolders
files    = dir();
names    = {files.name};
% Get a logical vector that tells which is a directory.
dir_flags = [files.isdir] & ~strcmp(names, '.') & ~strcmp(names, '..');
% Extract only those that are directories.
subdirs_names =  names(dir_flags) ;
rho_count = 0;
for subdir_num = 1:length( subdirs_names )
    subdir = char( subdirs_names( subdir_num ) )
    if strcmp(subdir, kr.input_dir)
        domain_loc = dir([subdir '/' kr.input_geom '*dat']);
        domain = load( [ subdir '/' domain_loc.name ] );
        domain = reshape(domain, flip(kr.domain_size));
        domain = erase_voxels(domain,kr.mesh_added,kr.num_slices);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        geom.name       = kr.input_geom;
        geom.print_size = true;
        geom.add_mesh   = false;
        geom.num_slices = 1; % add an additional slice to place the bc
        geom.swapXZ     = true;
        geom.scale_2    = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        create_geom_edist(logical(domain)*1,geom);
        
        vol_total = sum(~domain(:)); %total volume available for flow
        porosity = vol_total/length(domain(:));
        fprintf(fileID,'The sample porosity is %f \n',porosity);
    else
        fprintf(fileID,'---------------------- \n');
        fprintf(fileID,'Directory %s \n', subdir);
        fprintf(fileID,'---------------------- \n');
        rho_loc = dir([subdir '/rho*dat']);
        for rhos = 1:length(rho_loc)
            rho_count = rho_count + 1;
            rho = load( [ subdir '/' rho_loc.name] );
            rho = reshape(rho, flip(kr.domain_size));
            rho = erase_voxels(rho,kr.mesh_added,kr.num_slices);
            fluid_mask = rho>1;
            saturation = sum(fluid_mask(:))/vol_total;
            breakthrough = sum(sum( rho(:,:,end) )) > 1;
            
            %vel = load( [ in_folder all_vels(pc_inc).name] );
            %vel = reshape( vel, [3 flip(domain_size)] );
            %vz = squeeze( vel(1,:,:,:) ); %grab the first component of velocity (z)
            %%clear vel
            %vz = erase_voxels(vz,mesh_added,num_slices);
            %mean_vz  = mean( vz(fluid_mask)  )
            %mean_rho = mean( rho(fluid_mask) )
            %ca       = 0.1667*mean_vz/( 0.2 )
            
            % final domain for relative permeability
            domain4kr_f1 = ~fluid_mask;
            domain4kr_f2 = logical(domain + fluid_mask)*1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            name_kr         = ['_for_kr_' num2str(rho_count)];
            geom.print_size = false;
            geom.add_mesh   = false;
            geom.num_slices = 1; % add an additional slice to place the bc
            geom.swapXZ     = true;
            geom.scale_2    = false;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            geom.name = ['f1' name_kr];
            palabos_kr_f1 = create_geom_edist(domain4kr_f1,geom);
            
            geom.name = ['f2' name_kr];
            palabos_kr_f2 = create_geom_edist(domain4kr_f2,geom);
            
            fprintf(fileID,'--------------------------- \n');
            fprintf(fileID,'Simulation num %d \n',rho_count);
            
            if breakthrough == true
                fprintf(fileID,'Breakthrough has occured \n');
            else
                fprintf(fileID,'Breakthrough has NOT occured \n');
            end
            fprintf(fileID,'The saturation of f1 is %f \n',saturation);
        end
    end
    
end

fclose(fileID);

function array = erase_voxels(array,mesh_added,num_slices)
array(:,:,1:num_slices                       )=[];  %removes the first slices
array(:,:,end-num_slices-(mesh_added*2)+1:end)=[];  %removes the last slices + mesh +1 more
end

