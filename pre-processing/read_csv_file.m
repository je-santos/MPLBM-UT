function [imageVol]=read_csv_file(address, res)
%% convert sphere pack with center coordinates and radii to geometry for LBM simulations

% read csv file
files = dir(address);
infile1=files(1).name;
foldername= extractBefore(address,'*.');
data = readmatrix([foldername infile1]);
x_all = data(:,1);  % x-coordinate
y_all = data(:,2);  % y-coordinate
z_all = data(:,3);  % z-coordinate
r_all = data(:,4);  % radius
number=length(x_all);
max_X=max(x_all);

% Create cubical domain
MaxL=round(max_X*res);   %multiply by resolution for scaling
simSpace=zeros(MaxL,MaxL,MaxL);
simSizeX=MaxL;
simSizeY=MaxL;
simSizeZ=MaxL;
simSpace=uint8(simSpace);

% Adding spheres
for i=1:number
    
    r=res.*r_all(i);  %multiply by resolution for scaling 
    x=res.*x_all(i);
    y=res.*y_all(i);
    z=res.*z_all(i);
    
     % make 3d object
     [X,Y,Z]=meshgrid(-r:r,-r:r,-r:r);
     OBJ=(X/r).^2 + (Y/r).^2 + (Z/r).^2;
     OBJ=uint8(OBJ<1);
       
    [lx,ly,lz]=size(OBJ);
  
    x1=floor(x-lx/2)+1;
    x2=x1+lx-1;
    y1=floor(y-ly/2)+1;
    y2=y1+ly-1;
    z1=floor(z-lz/2)+1;
    z2=z1+lz-1;
    
    %crop if outside cubical domain
    
    if x1<1
        cutx=-(x1-1);
        x1=1;
        OBJ=OBJ(cutx+1:lx,:,:);
    end
    if x2>simSizeX
        cutx=x2-simSizeX;
        x2=simSizeX;
        OBJ=OBJ(1:lx-cutx,:,:);
    end
    if y1<1
        cuty=-(y1-1);
        y1=1;
        OBJ=OBJ(:,cuty+1:ly,:);
    end
    if y2>simSizeY
        cuty=y2-simSizeY;
        y2=simSizeY;
        OBJ=OBJ(:,1:ly-cuty,:);
    end
    if z1<1
        cutz=-(z1-1);
        z1=1;
        OBJ=OBJ(:,:,cutz+1:lz);
    end
    if z2>simSizeZ
        cutz=z2-simSizeZ;
        z2=simSizeZ;
        OBJ=OBJ(:,:,1:lz-cutz);
    end
    
     simSpace(x1:x2,y1:y2,z1:z2)=(simSpace(x1:x2,y1:y2,z1:z2)+OBJ);         
end

imageVol=uint8((simSpace)>0); %Convert to binary



 