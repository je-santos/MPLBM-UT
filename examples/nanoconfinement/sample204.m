% convert sphere pack with center coordinates and radii to image slices
% Abhishek Bihani, March 2019
% modified by Javier E. Santos, 2021
% Input Data can be found here: https://www.digitalrocksportal.org/projects/204


res = (0.6*10^2)*4; % For Scaling
ratios  = [10,7,5,2];
vol_frac = [30, 50, 49.9, 70];

for grain_r = ratios
    for vf = vol_frac
        data_files = dir(['raw_volumes\204_00\' num2str(grain_r) ...
                            '*_VL' num2str(vf) '*_1.csv']);
                        
        if numel(data_files)==0
            continue
        end
        
        for ind_file=1:numel(data_files)
            data = load(['raw_volumes\204_00\' data_files(ind_file).name]);
            csv_name = data_files(ind_file).name(1:end-4);
            if numel(csv_name)>12
                fix_str = csv_name(8:11);
                new_str = ceil(str2num(fix_str));
                csv_name = strrep(csv_name, fix_str, num2str(new_str));
            end
            
            x_all = data(:,1);  % x-coordinate
            y_all = data(:,2);  % y-coordinate
            z_all = data(:,3);  % z-coordinate
            r_all = data(:,4); % radius
            number=length(x_all);
            
            
            a = min(r_all);
            b = max(r_all);
            a1 = nnz(r_all==a);
            b1 = nnz(r_all==b);
            Vl= b1*4/3*pi*b^3;
            Vs= a1*4/3*pi*a^3;
            VL= Vl*100/(Vs+Vl);
            VL = round(VL,2);  % Volume percent of large grains
            ratio = b/a;       % Ratio of radius of large to small grains
            ratio=round(ratio);
            
            max_X=max(x_all);
            
            filename = sprintf('sliced_grainpack_%gto1_VL%g_',ratio,VL);% create output name
            
            %%%%%%%%
            
            MaxL=round(max_X*res);
            simSpace=zeros(MaxL,MaxL,MaxL);
            simSizeX=MaxL;
            simSizeY=MaxL;
            simSizeZ=MaxL;
            simSpace=uint8(simSpace);
            
            
            %Adding spheres
            for i=1:number
                
                r=(res.*r_all(i));
                x=res.*x_all(i);
                y=res.*y_all(i);
                z=res.*z_all(i);
                
                %make 3d object
                [X,Y,Z]=meshgrid(-r:r,-r:r,-r:r);
                OBJ=(X/r).^2 + (Y/r).^2 + (Z/r).^2;
                % OBJ=((X-xc)/rx).^2 + ((Y-yc)/ry).^2 + ((Z-zc)/rz).^2;
                OBJ=uint8(OBJ<1);
                
                rotObj=OBJ;
                
                %[rotObj,unused] = rotateDiscObj(OBJ,xRot,yRot,zRot);
                
                [lx,ly,lz]=size(rotObj);
                
                x1=floor(x-lx/2)+1;
                x2=x1+lx-1;
                y1=floor(y-ly/2)+1;
                y2=y1+ly-1;
                z1=floor(z-lz/2)+1;
                z2=z1+lz-1;
                %crop if outside sim grid
                
                if x1<1
                    cutx=-(x1-1);
                    x1=1;
                    rotObj=rotObj(cutx+1:lx,:,:);
                end
                if x2>simSizeX
                    cutx=x2-simSizeX;
                    x2=simSizeX;
                    rotObj=rotObj(1:lx-cutx,:,:);
                end
                if y1<1
                    cuty=-(y1-1);
                    y1=1;
                    rotObj=rotObj(:,cuty+1:ly,:);
                end
                if y2>simSizeY
                    cuty=y2-simSizeY;
                    y2=simSizeY;
                    rotObj=rotObj(:,1:ly-cuty,:);
                end
                if z1<1
                    cutz=-(z1-1);
                    z1=1;
                    rotObj=rotObj(:,:,cutz+1:lz);
                end
                if z2>simSizeZ
                    cutz=z2-simSizeZ;
                    z2=simSizeZ;
                    rotObj=rotObj(:,:,1:lz-cutz);
                end
                
                simSpace(x1:x2,y1:y2,z1:z2)=(simSpace(x1:x2,y1:y2,z1:z2)+rotObj.*10);
            end
            simSpace=uint8((simSpace)>0); %Convert to binary
            
            
            simSpace = simSpace( (1200-240):(1200+240-1), ...
                                 (1200-240):(1200+240-1), ...
                                 (1200-240):(1200+240-1) );
            
            save(['raw_volumes\204_00\voxelized\' csv_name],'simSpace')
            figure;imagesc(simSpace(:,:,240));
            pause(.1)
        end
    end
end




