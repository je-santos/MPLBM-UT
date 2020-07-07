addpath ('../../post-processing') %post-procesing libraries

%% Input for the function
kr.domain_size  = [200,200,200];
kr.mesh_added   = false; % was a neutral-wet mesh added at the end of the domain?
kr.num_slices   = 2;    % how many n empty slices at the begining and end of domain
kr.input_dir    = 'input';
kr.input_geom   = 'spheres4Palabos';
kr.output_file  = 'tmp/perc_path_info.txt';

percolation_path_calc