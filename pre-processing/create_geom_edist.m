function [geom4palabos]=create_geom_edist(data,name,num_slices,add_mesh,swapXZ,scale_2)

%%%% Inputs:
% data: image, where the pore-space is represented with zeros
% name: string with the name of the file for printing

tic

% Swap x and z data if needed to ensure Palabos simulation in Z-direction
if swapXZ == true
 data=permute(data,[3 2 1]);
end


%creates the computationally efficent geometry
edist = bwdist(~data);
geom4palabos = edist;
geom4palabos([1,end],:,:)=1; %makes sure that the outer box boundaries have
geom4palabos(:,[1,end],:)=1; %a wetting BC to avoid problems
geom4palabos(:,:,[1,end])=1;
geom4palabos(edist==0)=0;

geom4palabos(geom4palabos>0 & geom4palabos<2)=1;
geom4palabos(geom4palabos>1)=2;

% Double the grain (pore) size if needed to prevent single pixel throats
% for tight/ low porosity geometries
if scale_2 == true
    geom4palabos = imresize3(geom4palabos, 2, 'nearest'); 
end

% add a mesh if requested
if add_mesh == true
    mesh = toeplitz(mod(1:size(geom4palabos,2),2),mod(1:size(geom4palabos,3),2));
    mesh(mesh==1)=4; %change the mesh to be neutral-wet
    mesh_s = zeros([1 size(mesh)]);
    mesh_s(mesh==4)=4;
    geom4palabos = cat(1, geom4palabos, mesh_s);
end

% add a additional slices if requested
blank_slice = squeeze( geom4palabos(1:num_slices,:,:)*0 );
geom4palabos = cat(1, blank_slice, geom4palabos);
geom4palabos = cat(1, geom4palabos, blank_slice);

% print the geometry 
geom_name = [name '_' num2str( size(geom4palabos,1) ) '_' ...
     num2str( size(geom4palabos,2) ) '_' num2str( size(geom4palabos,3) ) '.dat'];
 
fid = fopen(['input/' geom_name], 'w'); % open the output file to write in

for x_coor=1:size(geom4palabos,1)
    fprintf(fid, '%i\n', squeeze( geom4palabos(x_coor,:,:) ) );
end

fclose(fid);

toc
end