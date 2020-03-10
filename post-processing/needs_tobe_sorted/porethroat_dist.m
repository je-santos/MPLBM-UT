        %Reading density vtk files generated from Palabos
        list=dir('vtkgeometry1.vti');
        listcell=struct2cell(list);
        image_name = listcell{1};
        
        image_ts(1)=int16(str2double(image_name));
        f.f1_vti_struct=xml2struct(image_name); %read output file
        %cd ../..
        f.f1_vti_str=base64decode(f.f1_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
        % f.f2_vti_str=base64decode(f.f2_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
        f.f1_vti_no=typecast([0 0 f.f1_vti_str],'single');
        % f.f2_vti_no=typecast([0 0 f.f2_vti_str],'double');
        vti_size=str2num(f.f1_vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
        vti_x=vti_size(1)+vti_size(2)+1;
        vti_y=vti_size(3)+vti_size(4)+1;
        vti_z=vti_size(5)+vti_size(6)+1;
        rhof1=reshape(f.f1_vti_no(3:end),[vti_x vti_y vti_z]);
        clear f
        tic
        
        
        % rhof1=simSpace;
        % rhof1=1-rhof1;  
        % k1=rhof1(:,:,60);
		% imagesc(k1);

		%Specify how many slices to remove from start and end along each axis
		remove = 125;
        
        rhof1(:,1:remove,:)=[];
        rhof1(:,end-remove:end,:)=[];
        rhof1(:,:,1:remove)=[];
        rhof1(:,:,end-remove:end)=[];
        rhof1(1:remove,:,:)=[];
        rhof1(end-remove:end,:,:)=[];
        
		rhof1(rhof1==0)=1;
        rhof1(rhof1<1)=0;
        
        porespace=rhof1(:,:,60);
        figure();
        imagesc(porespace);
        porespace1=squeeze(rhof1(:,60,:));
        figure();
        imagesc(porespace1);
		
        skel=Test_Skeleton3D(rhof1);       
        [NODE,LINK,perc_length] = Test_Skel2Graph3D(skel,rhof1);
        save('Pores','NODE');
        save('Throats','LINK');
        toc