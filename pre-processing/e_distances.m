function [e_domain,e_full, e_z] = e_distances( domain, o_print, num )

d_size = size(domain,1);
tmp = domain;
e_domain = bwdist(tmp);

for i = 1:d_size
    e_z(:,:,i) = bwdist( tmp(:,:,i) ); %Euclidean distance in Z
end

e_solid  =  bwdist(~tmp); %calculates the e dist in the solid
e_solid = -1*e_solid;
e_full = e_domain;
e_full(e_full==0) = e_solid(e_full==0);
clear tmp;

if o_print==true
    global saveto
    e_domain   = single(e_domain);
    e_z        = single(e_z);
    e_full     = single(e_full);
    save([ saveto '/euclidean_pore_',  num2str(num),'.mat'],'e_domain','-v7.3','-nocompression')
    save([ saveto '/euclidean_poreZ_', num2str(num),'.mat'],'e_z','-v7.3','-nocompression')
    save([ saveto '/euclidean_total_', num2str(num),'.mat'],'e_full','-v7.3','-nocompression')
end

end

