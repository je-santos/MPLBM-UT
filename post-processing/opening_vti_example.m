% Script to decode a .vti (palabos output) into a matlab matrix
% I gratefully acknowledge the help of Guillaume Flandin 

vti_struct=xml2struct('rho_f1_001_00291000.vti'); %read output file
vti_str=base64decode(vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
vti_no=typecast([0 0 vti_str],'single');
vti_size=str2num(vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
vti_x=vti_size(1)+vti_size(2)+1;
vti_y=vti_size(3)+vti_size(4)+1;
vti_z=vti_size(5)+vti_size(6)+1;


%vti_matrix=reshape(vti_no(2:end),[3 vti_x vti_y vti_z]);  %for velocity
vti_matrix=reshape(vti_no(2:end),[vti_x vti_y vti_z]);  %for density