addpath ('nano_lbm\') %pre-precesing libraries

sample_type = 'syn';

folder_loc =  ['matlab_volumes_' sample_type];
images = dir([folder_loc '/*.mat']);

short_name = 'tmp';
for i=1:numel(images) 
   name =  images(i).name(1:end-4);
   %disp(['374_' name])
   %disp([name])
   
   %if strcmp(short_name, name(1:end-4)) == false
   %    short_name = name(1:end-4);
   %    disp(['374_' short_name])
   %end
   
   im = unpackStruct(load([images(i).folder '/' images(i).name]));
   if isa(im, 'single')
      im = uint8(im); 
   end
   bcs = unique(im(:,:,1))+unique(im(:,:,end));
   if bcs>0
       disp([name ' does not percolate :('])
       a=1;
   end
   
   
end
