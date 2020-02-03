%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Hasan Khan
% UT-PGE 2016
% Code to create spherepack 
% with vuggy inclusions (optional)
% from XYZ coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all
tic 

%% options
write_dat_file          = false;
remove_boundary_spheres = false;


%% Input image
filename = 'FinneyPack.txt';
%filename = 'PeriodicPack_Matt_Balhoff_100.txt';
data = importdata(filename);
sphnum = data(:,1);
x = data(:,2); %[-18.38,18.31]
y = data(:,3); %[-17.97,17.71]
z = data(:,4); %[-18.24,18.16]
r = 1;
fprintf('Data loading complete\n')

%% Boundary size and vug removal

boundarysize = 7;
vugsquare = 0;

xmin = -vugsquare;
xmax =  vugsquare;
ymin = -vugsquare;
ymax =  vugsquare;
zmin = -vugsquare;
zmax =  vugsquare;
del = [];

for i = 1:numel(x)
    % ranges between [-boundarysize,boundarysize]
    if x(i)<-boundarysize || x(i)>boundarysize || y(i)<-boundarysize || ...
            y(i)>boundarysize || z(i)<-boundarysize || z(i)>boundarysize
        del = [del;i]; %#ok<*AGROW>
    end
    if x(i)>xmin && x(i)<xmax && y(i)>ymin && y(i)<ymax && ...
        z(i)>zmin && z(i)<zmax
        del = [del;i]; 
    end
end

if remove_boundary_spheres == true
    x(del) = []; y(del) = []; z(del) = []; % remove out of range spheres
    x = x + boundarysize; % shift the centers to comply  
    y = y + boundarysize; % with new range 
    z = z + boundarysize; % [0,2*boundarysize]
    fprintf('Vug zone created\n')
end

%% Generate sphere with specific radius and generate binary figure
scaling = 25;
rad = r * scaling;
SE = strel('sphere',rad);
sphere = SE.Neighborhood;

boxdim = 2*(boundarysize-r/2)*rad;
img = zeros(boxdim,boxdim,boxdim) > 0;
fprintf('Sphere neighbourhood generated\n')

%% Seed points for sphere centers

for i = 1:numel(x)
    xnew = ceil(x(i)*scaling);
    ynew = ceil(y(i)*scaling);
    znew = ceil(z(i)*scaling);
    img(xnew,ynew,znew) = 1;
    fprintf("Seed point generated # %d/%d\n",i,numel(x))
end

dilatedBW = imdilate(img,SE); % incorporate the spheres in the matrix
fprintf('\n%d spheres incorporated\n',numel(x))
matrixmesh = 1 - dilatedBW;

%% Output binary image

cropsize1 = floor(scaling/2);
cropsize2 = ceil(boxdim+scaling/2);
outputdims = cropsize2-cropsize1;
output = matrixmesh(cropsize1:cropsize2-1,cropsize1:cropsize2-1,...
    cropsize1:cropsize2-1);
fprintf('Writing image to RAW files\n')
fname = ['output_nx',num2str(outputdims),'_ny',num2str(outputdims),'_nz',...
    num2str(outputdims),'_littleEndian_8bit.raw'];
% fid=fopen(fname,'w+');
% cnt=fwrite(fid,output,'uint8','ieee-le');
% fclose(fid);

fprintf('(%i x %i x %i) image written with vug size = %i \n',...
    outputdims,outputdims,outputdims,vugsquare)

%% STL output

bin_out = output > 0.5;

%fv = isosurface(~output, 0.5);
%fprintf('Writing STL mesh file\n')
%stlwrite('mesh.stl',fv)
toc