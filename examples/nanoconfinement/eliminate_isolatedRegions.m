function [M,phi_f] = eliminate_isolatedRegions( M, connect )
% M is a binary image. Where 0 represents the pore space

tmp   = bwconncomp( ~M ,connect ); %check for unconnected regions in the pore-space
phi_i = 1-sum(M(:))/numel(M);  %Initial porosity


size_regions = cellfun(@numel,tmp.PixelIdxList);
[~, position] = sort(size_regions, 'descend');

if tmp.NumObjects>0
    for region=2:tmp.NumObjects
        M(tmp.PixelIdxList{position(region)})=1; %fill with 1's the unconnected regions
    end
end

phi_f    = 1-sum(M(:))/numel(M);  %porosity after conecting regions
phi_diff = phi_i-phi_f;


disp(['The final porosity is ' num2str(phi_f) ' where ' num2str(phi_diff) ...
      ' was erased. ' num2str(tmp.NumObjects) ' regions were found.'])
end

