function [tOf_L,tOf_R] = tOf( domain, o_print, num )

d_size      = size(domain,1);
SpeedImage  = double(domain);
X           = repmat(1:d_size,1,d_size);
Y = [];

for i=1:d_size
    Y(end+1:end+d_size) = i; %injection points
end

SourcePoint_L = [ X ; Y ; ones(1,numel(X))        ]; %inject from the Left
SourcePoint_R = [ X ; Y ; ones(1,numel(X))*d_size ]; %inject from Right

SpeedImage( domain==0 ) = 1e-4; %give the solids a low conductivity

tic; tOf_L = msfm( SpeedImage, SourcePoint_L, true, true ); toc;
%tic; tOf_R = msfm( SpeedImage, SourcePoint_R, true, true ); toc;

tmp     = linspace( 0, d_size-1, d_size);
linear  = zeros(d_size, d_size, d_size );
for i=1:d_size
    linear(:,:,i) = ones( d_size, d_size )*tmp(i); %computes dist from inlet-1
end
clear tmp

% detrend
tOf_L = tOf_L - linear;
%tOf_R = tOf_R - flip(linear,3);

tOf_L( domain==0 ) = 0; %erase the solids
%tOf_R( domain_inv==0 ) = 0; %erase the solids

max_time = 5000; %check histogram first
fprintf('L: The num of vals>threshold is %d \n',sum(sum( sum( tOf_L >  max_time ) )))
%fprintf('R: The num of vals>threshold is %d \n',sum(sum( sum( tOf_R >  max_time ) )))
fprintf('L: The num of vals<threshold is %d \n',sum(sum( sum( tOf_L < -max_time ) )))
%fprintf('R: The num of vals<threshold is %d \n',sum(sum( sum( tOf_R < -max_time ) )))

tOf_L( tOf_L > max_time ) = 0;
%tOf_R( tOf_R > max_time ) = 0;

%tOf_L( tOf_L < -max_time ) = 0;
%tOf_R( tOf_R < -max_time ) = 0;


figure;
subplot(1,2,1)
imagesc(squeeze( domain( round(d_size/2) ,:,:) ));
colorbar
title('Fluid velocity')
subplot(1,2,2)
imagesc(squeeze( tOf_L( round(d_size/2) ,:,:) ));
colorbar
title('Time of Flight')


if o_print==true
    global saveto
    tOf_L = single(tOf_L);
    %tOf_R = single(tOf_R);
    save([ saveto '/ToF_l_' , num2str(num),'.mat'],'tOf_L','-v7.3','-nocompression')
    %save([ saveto '/ToF_r_' , num2str(num),'.mat'],'tOf_R','-v7.3','-nocompression')
end

end

