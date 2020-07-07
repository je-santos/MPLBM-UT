addpath ('../../pre-processing') %pre-procesing libraries

folder = fileparts(which(mfilename)); % Determine location of file folder (post-processing)
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geom.name       = kr.input_geom;
geom.print_size = true;
geom.add_mesh   = false;
geom.num_slices = 0; % add an additional slice to place the bc
geom.swapXZ     = true;
geom.scale_2    = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
first = 0;
for subdir_num = 1:length( subdirs_names )
    subdir = char( subdirs_names( subdir_num ) )
  %  if ~strcmp(subdir, kr.input_dir)
        fprintf(fileID,'---------------------- \n');
        fprintf(fileID,'Directory %s \n', subdir);
        fprintf(fileID,'---------------------- \n');
        rho_loc = dir([subdir '/rho*dat']);
        
        for rhos = 1:length(rho_loc)
            rho_count = rho_count + 1;
            rho = load( [ subdir '/' rho_loc(rhos).name] );
            rho = reshape(rho, kr.domain_size);
            %rho = erase_voxels(rho,kr.mesh_added,kr.num_slices);
            Front_f1 = 0;
            
            if geom.swapXZ  == true
                rho = permute(rho,[3 2 1]);
            end
            if geom.scale_2 == true
                rho = imresize3(rho, 2, 'nearest');
            end
            
            rhof1 = rho;
            rhof1(rhof1==0.4 | rhof1==-0.4)=0;  %Converting boundary pixels to grains
            rho1_temp=rhof1;
            rho1_temp(rho1_temp>1)=1;
            rho1_temp(rho1_temp<1)=0;
            
            %Find movement of first non-wetting front/finger percolating in X-direction
            rho1_temp2=rho1_temp;
            rho1_temp2(1:kr.num_slices,:,:)=[];
            rho1_temp2(end-kr.num_slices:end,:,:)=[];
            X=size(rho1_temp2,1);
            
            saturation = sum(sum(sum(rhof1>1)))/sum(sum(sum(rhof1~=0)));
            
            
            for count=1:X
                rho1_temp3=squeeze(rho1_temp2(count,:,:));
                leng_f1(count) = max(max(rho1_temp3));
            end
            idx = find(leng_f1 == 1);
            Front_f1=max(idx);
                       
            fprintf(fileID,'--------------------------- \n');
            fprintf(fileID,'Simulation num %d \n',rho_count);
            
            if Front_f1 == X && first == 0 % Non-wetting fluid has percolated
                fprintf(fileID,'Breakthrough has occured for simulation num %d \n',rho_count);
                
                %Find the skeleton of the invading fluid at breakthrough and calculates percolation path.
                skel=Test_Skeleton3D(rho1_temp);
                [NODE,LINK,perc_length] = Test_Skel2Graph3D(skel,rho1_temp);
                
                perc_length = round(perc_length,3);
                
                tortuosity = perc_length/X;
                tortuosity = round(tortuosity,3);
                
                first=1;
                
                fprintf(fileID,'Percolation path length %d and domain length in flow direction %d \n',perc_length,X);
                fprintf(fileID,'Percolation path tortuosity at breakthrough is %d \n',tortuosity);
                
            end
            fprintf(fileID,'The saturation of f1 is %f \n',saturation);
        end
  %  end   
end

fclose(fileID);
