%This code will calculate the wetting saturation for all vtk files, will convert
%the fluid configurations (1 and 2) to .dat files for 1-phase LBM
%simulation and will find the vtk file where breakthrough occurs and the
%percolation path/tortuosity. You can choose if you want to generate fluid
%geometries or calculate percolation path at breakthrough

clear
tic
%Options:

medial_axis=0; %To find fluid medial axes and first percolation path
print=0;      %Create geometries for fluids 1 and 2 for relative permeability calculations


%post-processing new...

%Update path to folder where output .dat files from 2-phase simulation are stored
directory='D:\LBM C\MPLBM-UT\post-processing\input_files\';

list=dir([directory 'rho_f1*.dat']);
listcell=struct2cell(list);
[a,b]=size(listcell);
listcell=listcell(1,:);
first=0;

%Reading 1 vtk file generated from Palabos to get geometry limits
f.f1_vti_struct=xml2struct([directory 'porousMedium.vti']); %read output file
f.f1_vti_str=base64decode(f.f1_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
f.f1_vti_no=typecast([0 0 f.f1_vti_str],'double');
vti_size=str2num(f.f1_vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
z_lim=vti_size(1)+vti_size(2)+1;
y_lim=vti_size(3)+vti_size(4)+1;
x_lim=vti_size(5)+vti_size(6)+1;

% Looping through capillary pressure increments to calculate Sw and to
% generate geometries.
for I=1:b
     
    %Reading density .dat files generated from Palabos
    image_name = listcell{I};
    image_ts=load([directory image_name]);
    rhof1 = reshape(image_ts, x_lim, y_lim, z_lim);
      
    
    %Inverting X and Z directions to keep geometry consistent
    rhof1 = permute(rhof1,[3 2 1]);
    
    
    % Identifying only Fluid 1 and Fluid 2 in entire geometry
    rhof1(rhof1==0.4 | rhof1==-0.4)=0;  %Converting boundary pixels to grains
    rho1_temp=rhof1;
    rho1_temp(rho1_temp>1)=1;
    rho1_temp(rho1_temp<1)=0;
    
    rho2_temp=rhof1;
    rho2_temp(rho2_temp>1)=0;
    rho2_temp(rho2_temp<1 & rho2_temp>0)=1;
    
    
    %Find movement of first non-wetting front/finger percolating in X-direction
    rho1_temp2=rho1_temp;
    rho1_temp2(1:2,:,:)=[];
    rho1_temp2(end-2:end,:,:)=[];
    X=size(rho1_temp2,1);
    
    for count=1:X
        rho1_temp1=squeeze(rho1_temp2(count,:,:));
        leng_f1(count) = max(max(rho1_temp1));
    end
    idx = find(leng_f1 == 1);
    Front_f1(I)=max(idx);
    
    %Find percolation length and tortuosity of first breakthrough path
    
    if medial_axis==1
        if Front_f1(I) == X && first == 0 % Non-wetting fluid has percolated
            skel=Test_Skeleton3D(rho1_temp);
            [NODE,LINK,perc_length] = Test_Skel2Graph3D(skel,rho1_temp); %This function 
            % finds the skeleton of the invading fluid at breakthrough and calculates percolation path.
            tortuosity = perc_length/X;
            first=1;
            disp('Breakthrough occurs at: ');
            disp(image_name);
            Breakthrough=image_name;
            breakthrough=[tortuosity,perc_length];
            csvwrite([directory '\' 'breakthrough' '.csv'],breakthrough);
        end
    end
    
    %Create geometries for fluids 1 and 2 for relative permeability calculations
    
    
    if print==1
        
        for count=1:2
            if count==1
               data = rho1_temp;
            %    data = permute(rho1_temp,[3 2 1]);
                baseOutput = sprintf('lattice_f1_forK_%g_',I);% create output name
            else
               data = rho2_temp;
            %    data = permute(rho2_temp,[3 2 1]);
                baseOutput = sprintf('lattice_f2_forK_%g_',I);% create output name
            end
            
            numFiles=size(data,3);
            
            baseOutput = [directory baseOutput '.dat'];
            
            
            %creates the computationally efficent geometry
            edist = bwdist(data);
            geom4palabos = edist;
            geom4palabos([1,end],:,:)=1; %makes sure that the outer box boundaries have
            geom4palabos(:,[1,end],:)=1; %a wetting BC to avoid problems
            geom4palabos(:,:,[1,end])=1;
            geom4palabos(edist==0)=0;
            
            geom4palabos(geom4palabos>0 & geom4palabos<2)=1;
            geom4palabos(geom4palabos>1)=2;
                       
            fid = fopen(baseOutput, 'w');    % open the output file to write in
            
            for x_coor=1:size(geom4palabos,1)
                fprintf(fid, '%i\n', squeeze( geom4palabos(x_coor,:,:) ) );
            end
            
            fclose(fid);
            
        end
    end
    R=1;
    sat_nw(R,I)=sum(sum(sum(rhof1>1)))/sum(sum(sum(rhof1~=0)));
    sat_w(R,I)= 1 - sat_nw(R,I);
    csvwrite([directory '\' 'Sw' '.csv'],sat_w); % Saves wetting saturation for each pressure increment
    clear rhof1
    clear rho1_temp
    clear rho2_temp
end
toc


