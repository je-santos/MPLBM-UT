function [ M ] = eliminate_isolatedRegions( M, connect )
% M is a binary image. Where 0 represents the pore space and 1 the solid

tmp=bwconncomp( ~M ,connect ); %check for unconnected regions in the pore-space
phi_i=1-sum(M(:))/numel(M);  %Initial porosity

fprintf('The intial porosity is %.3f %% with %d connected regions \n', ...
            phi_i*100, tmp.NumObjects)


size_regions = cellfun(@numel,tmp.PixelIdxList);
[ind_size, position] = sort(size_regions, 'descend');

if tmp.NumObjects>0
    for region=2:tmp.NumObjects
        M(tmp.PixelIdxList{position(region)})=1; %fill with 1's the unconnected regions
    end
end

phi_f=1-sum(M(:))/500^3;  %porosity after conecting regions
phi_diff=(phi_i-phi_f)*100;

fprintf('The porosity difference with the fully-connected image is %.3f %% \n', phi_diff)

end

