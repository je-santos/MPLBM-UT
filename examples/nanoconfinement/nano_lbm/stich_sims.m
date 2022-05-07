sample_type = 'real';

folder_loc =  ['../matlab_volumes_' sample_type];
output_loc = [folder_loc  '_tiff'];
images = dir([folder_loc '/*.mat']);




% Get MFP
tiff_dir = output_loc;
%mfp_loc = ['../mfp_' sample_type];
mfp_loc = ['../tmp'];
%sim_geom_loc = ['../domains_' sample_type];
sim_geom_loc = ['../big_domains'];


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
       tmp = split(folders(f).name, '_');
       sim_size = str2num( tmp{end-1} );
       if ( sim_size == 256 | sim_size == 480 ) == false
          error('Size is wrong') 
       end
       output_tensors_fromSim(output_dir,folders(f).name,sim_size,save_to)
    end
end


