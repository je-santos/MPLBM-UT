addpath ('../../pre-processing') %pre-precesing libraries
global saveto        
saveto = 'tmp';

domain = load('tmp/carbonate4Palabos_252_250_250_vel.dat');
domain_size = [252,250,250];
domain = reshape(domain, [3 flip(domain_size)]);
domain = squeeze( domain(1,:,:,:) );

domain(:,:,1  )=[];
domain(:,:,end)=[]; %remove empty slices

min_vel   = min(  abs(domain(domain~=0)));
range_vel = range(abs(domain(domain~=0)));
domain(domain~=0) = (domain(domain~=0) - min_vel)/range_vel;

[tOf_L,tOf_R] = tOf( domain, true, 1 );
 
