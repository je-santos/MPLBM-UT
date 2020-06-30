%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to create a flat surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
dom.name='flat_surface';
dom.size_x=100;
dom.size_y=100;
dom.size_z=100;


fid_domain1=fopen([dom.name '.dat'],'w');
C=ones(dom.size_x,dom.size_y);  % floor
B=C*0;                          % empty space

for i=1:dom.size_z  %prints in Palabos depth
    if i==1
        fprintf(fid_domain1,'%i\n',C);
    else
        fprintf(fid_domain1,'%i\n',B);
    end
   
end


fclose(fid_domain1);
disp('Flat Surface Created')