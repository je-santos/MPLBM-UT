%This code will calculate the wetting saturation for all vtk files, will convert
%the fluid configurations (1 and 2) to .dat files for 1-phase LBM
%simulation and will find the vtk file where breakthrough occurs and the
%percolation path/tortuosity. You can choose if you want to generate fluid
%geometries or calculate percolation path at breakthrough

clear

%Options:

medial_axis=0; %To find fluid medial axes and first percolation path
print=1;      %Create geometries for fluids 1 and 2 for relative permeability calculations 

%Update path to folder where output vti files from 2-phase simulation are stored
directory='C:\Users\Abhishek\Desktop\MultiphasePorousMediaPalabos-master\src\2-phase_LBM\spheres_test\'; 

    
    
    list=dir([directory 'rho_f1*.vti']);
    listcell=struct2cell(list);
    [a,b]=size(listcell);
    listcell=listcell(1,:);
    first=0;
    
  
    for I=1:b
        
        %Reading density vtk files generated from Palabos
        
        image_name = listcell{I};
        image_ts(I)=int16(str2double([directory image_name]));
        f.f1_vti_struct=xml2struct([directory image_name]); %read output file
        %cd ../..
        f.f1_vti_str=base64decode(f.f1_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
        % f.f2_vti_str=base64decode(f.f2_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
        f.f1_vti_no=typecast([0 0 f.f1_vti_str],'double');
        % f.f2_vti_no=typecast([0 0 f.f2_vti_str],'double');
        vti_size=str2num(f.f1_vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
        vti_x=vti_size(1)+vti_size(2)+1;
        vti_y=vti_size(3)+vti_size(4)+1;
        vti_z=vti_size(5)+vti_size(6)+1;
        rhof1=reshape(f.f1_vti_no(2:end),[vti_x vti_y vti_z]);
        clear f
        
         rhof1(1:2,:,:)=[];
        
  
        rhof1(end-2:end,:,:)=[];
        
        rhof1(rhof1==0.4 | rhof1==-0.4)=0;


% Identifying only Fluid 1 and Fluid 2 in entire geometry

        rho1_temp=rhof1;
        rho1_temp(rho1_temp>1)=1;
        rho1_temp(rho1_temp<1)=0;
        rho2_temp=rhof1;
        rho2_temp(rho2_temp>1)=0;
        rho2_temp(rho2_temp<1 & rho2_temp>0)=1;

        
%         sizef=size(rhof1);
        
%                 for j=1:sizef(1)
%                     image(squeeze(rhof1(j,:,:)),'CDataMapping','scaled')
%                     sat_nw_front(I,j)=sum(sum(squeeze(rhof1(j,:,:))>1))/sum(sum(sum(rhof1(j,:,:)~=0)));
%                     title(['Sw=',num2str(sat_nw_front(I,j))])
%                     colorbar
%                     pause(.2)
%                 end
        
        
      %Find movement of first non-wetting front/finger percolating in X-direction
      
        X=size(rho1_temp,1);
        for count=1:X
        rho1_temp1=squeeze(rho1_temp(count,:,:));
        leng_f1(count) = max(max(rho1_temp1));
        end
        idx = find(leng_f1 == 1);
        Front_f1(I)=max(idx);
       
     %Find percolation length and tortuosity   
     
     
     if medial_axis==1   
         if Front_f1(I) == X && first == 0 % Non-wetting fluid has percolated
             skel=Test_Skeleton3D(rho1_temp);
             [NODE,LINK,perc_length] = Test_Skel2Graph3D(skel,rho1_temp);
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
                  %   data = rho1_temp;
                    data = permute(rho1_temp,[3 2 1]);
                    baseOutput = sprintf('lattice_f1_forK_%g_',I);% create output name
                else
                 %   data = rho2_temp;
                    data = permute(rho2_temp,[3 2 1]);
                    baseOutput = sprintf('lattice_f2_forK_%g_',I);% create output name
                end
                
                numFiles=size(data,3);
                
                
                
                baseOutput = [baseOutput '.dat'];
                
                fid = fopen(baseOutput, 'w');    % open the output file to write in
                
                
                
                %%%%%%%%%%%%%%%%%%%%% INLET SLICE %%%%%%%%%%%%%%%%%%%%%%%%
                % fname = [basename num2str(1, '%0.4i') '.tif'];
                % fnameu = [basename num2str(2, '%0.4i') '.tif'];
                
                BB = data(:,:,1);
                CC = data(:,:,2);
                
                nx=size(BB,2);
                ny=size(BB,1);
                
                % imagesc(BB);
                
                B=zeros(ny,nx);
                
                wholeGeom=zeros(ny,nx,2);
                
                
                wholeGeom(:,:,1)=BB;
                wholeGeom(:,:,2)=CC;
                
                indexMin=find(wholeGeom==0);
                indexMax=find(wholeGeom>0);
                
                wholeGeom(indexMin)=255;
                wholeGeom(indexMax)=0;
                
                
                %i
                rA = circshift(wholeGeom,[0,1, 0]);
                lA = circshift(wholeGeom,[0,-1, 0]);
                fA = circshift(wholeGeom,[1,0, 0]);
                bA = circshift(wholeGeom,[-1,0,0]);
                
                rfA = circshift(wholeGeom,[1,1,0]);
                rbA = circshift(wholeGeom,[-1,1,0]);
                lfA = circshift(wholeGeom,[1,-1,0]);
                lbA = circshift(wholeGeom,[-1,-1,0]);
                
                
                %i+1
                uA=circshift(wholeGeom,[0,0, 1]);
                
                urA = circshift(wholeGeom,[0,1, 1]);
                ulA = circshift(wholeGeom,[0,-1, 1]);
                ufA = circshift(wholeGeom,[1,0, 1]);
                ubA = circshift(wholeGeom,[-1,0,1]);
                
                urfA = circshift(wholeGeom,[1,1,1]);
                urbA = circshift(wholeGeom,[-1,1,1]);
                ulfA = circshift(wholeGeom,[1,-1,1]);
                ulbA = circshift(wholeGeom,[-1,-1,1]);
                
                for i=1:nx
                    for j=1:ny
                        
                        if (wholeGeom(j,i,1) == 255 && rA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && lA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && fA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && bA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,1) == 255 && rfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && rbA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && lfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && lbA(j,i)==0)
                            B(j,i) = 1;
                            
                            
                            
                        elseif (wholeGeom(j,i,1) == 255 && uA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,1) == 255 && urA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && ulA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && ufA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && ubA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,1) == 255 && urfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && urbA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && ulfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,1) == 255 && ulbA(j,i)==0)
                            B(j,i) = 1;
                            
                            
                        elseif (wholeGeom(j,i,1) == 0)
                            B(j,i) = 0;
                        else
                            B(j,i)=2;
                            
                        end
                        
                    end
                end
                
                % image(30*B)
                % axis equal
                % drawnow
                
                'printing first slice'
                fprintf(fid,'%i\n',B*0);
                fprintf(fid,'%i\n',B*0);
                fprintf(fid, '%i\n', B);
                
                
                %%%%%%%%%%%%%%%%%%%%% INTERNAL SLICES %%%%%%%%%%%%%%%%%%%%%%%%
                
                
                for ii=2:numFiles-1
                    ii
                    % fname = [basename num2str(ii, '%0.4i') '.tif'];
                    % fnamed =[basename num2str(ii-1, '%0.4i') '.tif'];
                    % fnameu = [basename num2str(ii+1, '%0.4i') '.tif'];
                    
                    AA = data(:,:,ii-1);
                    BB = data(:,:,ii);
                    CC = data(:,:,ii+1);
                    
                    % AA=imread(fnamed ,'tif');
                    % BB=imread(fname ,'tif');
                    % CC=imread(fnameu ,'tif');
                    
                    wholeGeom=zeros(ny,nx,3);
                    
                    wholeGeom(:,:,1)=AA;
                    wholeGeom(:,:,2)=BB;
                    wholeGeom(:,:,3)=CC;
                    
                    indexMin=find(wholeGeom==0);
                    indexMax=find(wholeGeom>0);
                    
                    wholeGeom(indexMin)=255;
                    wholeGeom(indexMax)=0;
                    
                    %i
                    rA = circshift(wholeGeom,[0,1, 0]);
                    lA = circshift(wholeGeom,[0,-1, 0]);
                    fA = circshift(wholeGeom,[1,0, 0]);
                    bA = circshift(wholeGeom,[-1,0,0]);
                    
                    rfA = circshift(wholeGeom,[1,1,0]);
                    rbA = circshift(wholeGeom,[-1,1,0]);
                    lfA = circshift(wholeGeom,[1,-1,0]);
                    lbA = circshift(wholeGeom,[-1,-1,0]);
                    
                    
                    %i-1
                    dA=circshift(wholeGeom,[0,0, -1]);
                    
                    drA = circshift(wholeGeom,[0,1, -1]);
                    dlA = circshift(wholeGeom,[0,-1, -1]);
                    dfA = circshift(wholeGeom,[1,0, -1]);
                    dbA = circshift(wholeGeom,[-1,0,-1]);
                    
                    drfA = circshift(wholeGeom,[1,1,-1]);
                    drbA = circshift(wholeGeom,[-1,1,-1]);
                    dlfA = circshift(wholeGeom,[1,-1,-1]);
                    dlbA = circshift(wholeGeom,[-1,-1,-1]);
                    
                    %i+1
                    uA=circshift(wholeGeom,[0,0, 1]);
                    
                    urA = circshift(wholeGeom,[0,1, 1]);
                    ulA = circshift(wholeGeom,[0,-1, 1]);
                    ufA = circshift(wholeGeom,[1,0, 1]);
                    ubA = circshift(wholeGeom,[-1,0,1]);
                    
                    urfA = circshift(wholeGeom,[1,1,1]);
                    urbA = circshift(wholeGeom,[-1,1,1]);
                    ulfA = circshift(wholeGeom,[1,-1,1]);
                    ulbA = circshift(wholeGeom,[-1,-1,1]);
                    
                    for i=1:nx
                        for j=1:ny
                            if (wholeGeom(j,i,2) == 255 && rA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && lA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && fA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && bA(j,i)==0)
                                B(j,i) = 1;
                                
                            elseif (wholeGeom(j,i,2) == 255 && rfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && rbA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && lfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && lbA(j,i)==0)
                                B(j,i) = 1;
                                
                                
                                
                            elseif (wholeGeom(j,i,2) == 255 && dA(j,i)==0)
                                B(j,i) = 1;
                                
                            elseif (wholeGeom(j,i,2) == 255 && drA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && dlA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && dfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && dbA(j,i)==0)
                                B(j,i) = 1;
                                
                            elseif (wholeGeom(j,i,2) == 255 && drfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && drbA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && dlfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && dlbA(j,i)==0)
                                B(j,i) = 1;
                                
                                
                            elseif (wholeGeom(j,i,2) == 255 && uA(j,i)==0)
                                B(j,i) = 1;
                                
                            elseif (wholeGeom(j,i,2) == 255 && urA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && ulA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && ufA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && ubA(j,i)==0)
                                B(j,i) = 1;
                                
                            elseif (wholeGeom(j,i,2) == 255 && urfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && urbA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && ulfA(j,i)==0)
                                B(j,i) = 1;
                            elseif (wholeGeom(j,i,2) == 255 && ulbA(j,i)==0)
                                B(j,i) = 1;
                                
                                
                            elseif (wholeGeom(j,i,2) == 0)
                                B(j,i) = 0;
                            else
                                B(j,i)=2;
                                
                            end
                            
                        end
                    end
                    
                    % image(30*B)
                    % axis equal
                    % drawnow
                    
                    fprintf(fid, '%i\n', B);
                end
                
              
                %%%%%%%%%%%%%%%%%%%%% OUTLET SLICES %%%%%%%%%%%%%%%%%%%%%%%%
                
                % fname = [basename num2str(numFiles, '%0.4i') '.tif'];
                % fnamed =[basename num2str(numFiles-1, '%0.4i') '.tif'];
                
                % AA=imread(fnamed ,'tif');
                % BB=imread(fname ,'tif');
                
                AA = data(:,:,numFiles-1);
                BB = data(:,:,numFiles);
                
                wholeGeom=zeros(ny,nx,2);
                
                wholeGeom(:,:,1)=AA;
                wholeGeom(:,:,2)=BB;
                
                indexMin=find(wholeGeom==0);
                indexMax=find(wholeGeom>0);
                
                wholeGeom(indexMin)=255;
                wholeGeom(indexMax)=0;
                
                %i
                rA = circshift(wholeGeom,[0,1, 0]);
                lA = circshift(wholeGeom,[0,-1, 0]);
                fA = circshift(wholeGeom,[1,0, 0]);
                bA = circshift(wholeGeom,[-1,0,0]);
                
                rfA = circshift(wholeGeom,[1,1,0]);
                rbA = circshift(wholeGeom,[-1,1,0]);
                lfA = circshift(wholeGeom,[1,-1,0]);
                lbA = circshift(wholeGeom,[-1,-1,0]);
                
                
                %i-1
                dA=circshift(wholeGeom,[0,0, -1]);
                
                drA = circshift(wholeGeom,[0,1, -1]);
                dlA = circshift(wholeGeom,[0,-1, -1]);
                dfA = circshift(wholeGeom,[1,0, -1]);
                dbA = circshift(wholeGeom,[-1,0,-1]);
                
                drfA = circshift(wholeGeom,[1,1,-1]);
                drbA = circshift(wholeGeom,[-1,1,-1]);
                dlfA = circshift(wholeGeom,[1,-1,-1]);
                dlbA = circshift(wholeGeom,[-1,-1,-1]);
                
                
                
                for i=1:nx
                    for j=1:ny
                        if (wholeGeom(j,i,2) == 255 && rA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && lA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && fA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && bA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,2) == 255 && rfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && rbA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && lfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && lbA(j,i)==0)
                            B(j,i) = 1;
                            
                            
                            
                        elseif (wholeGeom(j,i,2) == 255 && dA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,2) == 255 && drA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && dlA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && dfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && dbA(j,i)==0)
                            B(j,i) = 1;
                            
                        elseif (wholeGeom(j,i,2) == 255 && drfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && drbA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && dlfA(j,i)==0)
                            B(j,i) = 1;
                        elseif (wholeGeom(j,i,2) == 255 && dlbA(j,i)==0)
                            B(j,i) = 1;
                            
                            
                            
                        elseif (wholeGeom(j,i,2) == 0)
                            B(j,i) = 0;
                        else
                            B(j,i)=2;
                            
                        end
                        
                    end
                end
         
                
                'printing last slice'
                fprintf(fid, '%i\n', B);
                
                [r,c]=size(B);
                

                
                fprintf(fid,'%i\n',B*0);
                fprintf(fid,'%i\n',B*0);
                fprintf(fid,'%i\n',B*0);
                fclose(fid);
                
            end
        end
        R=1;
        sat_nw(R,I)=sum(sum(sum(rhof1>1)))/sum(sum(sum(rhof1~=0)));
        sat_w(R,I)= 1 - sat_nw(R,I);
        csvwrite([directory '\' 'Sw' '.csv'],sat_w);
        clear rhof1
        clear rho1_temp
     %  save('rho1_temp.mat','rho1_temp')
        clear rho2_temp
    end



