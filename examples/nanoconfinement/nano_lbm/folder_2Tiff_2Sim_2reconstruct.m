sample_type = 'real';

folder_loc =  ['../matlab_volumes_' sample_type];
output_loc = [folder_loc  '_tiff'];
images = dir([folder_loc '/*.mat']);



for i=1:numel(images) 
   name =  images(i).name(1:end-4);
   im = unpackStruct(load([images(i).folder '/' images(i).name]));
   if isa(im, 'single')
      im = uint8(im); 
   end
   Vol2Tiff(im, output_loc, name)
end


% Get MFP
tiff_dir = output_loc;
%mfp_loc = ['../mfp_' sample_type];
mfp_loc = ['../tmp'];
%sim_geom_loc = ['../domains_' sample_type];
sim_geom_loc = ['../big_domains'];


parfor i=20:numel(images) 
   sample_name =  images(i).name(1:end-4)
   sample_size =  sample_name(end-2:end);
   domains_4sim(sample_name, tiff_dir, mfp_loc, sim_geom_loc);
end


%% Reconstruct LBM solution
output_dir = ['../domains_' sample_type];
save_to = ['../results_' sample_type];

pres = [1,2,5,10,20];
for p=pres
    folders = dir([output_dir '/*' num2str(p)]);
    dirFlags = [folders.isdir];
    folders = folders(dirFlags);
    
    for f=1:numel(folders)
       disp(folders(f).name)
       output_tensors_fromSim(output_dir,folders(f).name,256,save_to)
    end
end


