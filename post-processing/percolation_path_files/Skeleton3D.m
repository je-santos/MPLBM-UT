function skel = Skeleton3D(skel,spare)
% SKELETON3D Calculate the 3D skeleton of an arbitrary binary volume using parallel medial axis thinning.
%
% skel = SKELETON3D(img) returns the skeleton of the binary volume 'img'
% skel = SKELETON3D(img,mask) preserves foreground voxels in 'mask'
%
% MATLAB vectorized implementation of the algorithm by Lee, Kashyap and Chu
% "Building skeleton models via 3-D medial surface/axis thinning algorithms."
% Computer Vision, Graphics, and Image Processing, 56(6):462�478, 1994.
%
% Inspired by the ITK implementation by Hanno Homann
% http://hdl.handle.net/1926/1292
% and the Fiji/ImageJ plugin by Ignacio Arganda-Carreras
% http://fiji.sc/wiki/index.php/Skeletonize3D
%
% Philip Kollmannsberger (philipk@gmx.net)
%
% For more information, see <a
% href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/43400-skeleton3d')">Skeleton3D</a> at the MATLAB File Exchange.

% pad volume with zeros to avoid edge effects
% skel = logical(skel);
skel=padarray(skel,[1 1 1]);

if(nargin==2)
    spare=padarray(spare,[1 1 1]);
end

% fill lookup table
eulerLUT = FillEulerLUT;

width = size(skel,1);
height = size(skel,2);
depth = size(skel,3);

unchangedBorders = 0;

while( unchangedBorders < 6 )  % loop until no change for all six border types
    unchangedBorders = 0;
    for currentBorder=1:6 % loop over all 6 directions
        cands=false(width,height,depth, 'like', skel);
        switch currentBorder
            case 4
                x=2:size(skel,1); % identify border voxels as candidates
                cands(x,:,:)=skel(x,:,:) - skel(x-1,:,:);
            case 3
                x=1:size(skel,1)-1;
                cands(x,:,:)=skel(x,:,:) - skel(x+1,:,:);
            case 1
                y=2:size(skel,2);
                cands(:,y,:)=skel(:,y,:) - skel(:,y-1,:);
            case 2
                y=1:size(skel,2)-1;
                cands(:,y,:)=skel(:,y,:) - skel(:,y+1,:);
            case 6
                z=2:size(skel,3);
                cands(:,:,z)=skel(:,:,z) - skel(:,:,z-1);
            case 5
                z=1:size(skel,3)-1;
                cands(:,:,z)=skel(:,:,z) - skel(:,:,z+1);
        end
        
        % if excluded voxels were passed, remove them from candidates
        if(nargin==2)
            cands = cands & ~spare;
        end
        
        % make sure all candidates are indeed foreground voxels
        cands = cands(:)==1 & skel(:)==1;
        
        noChange = true;
                    
        if any(cands)
            cands = find(cands);
            % get subscript indices of candidates
            [x,y,z]=ind2sub([width height depth],cands);
            
            % get 26-neighbourhood of candidates in volume
            nhood = pk_get_nh(skel,cands);
            
            % remove all endpoints (exactly one nb) from list
            di1 = sum(nhood,2)==2;
            nhood(di1,:)=[];
            cands(di1)=[];
            x(di1)=[];
            y(di1)=[];
            z(di1)=[];
            
            % remove all non-Euler-invariant points from list
            di2 = ~p_EulerInv(nhood, eulerLUT);
            nhood(di2,:)=[];
            cands(di2)=[];
            x(di2)=[];
            y(di2)=[];
            z(di2)=[];
            
            % remove all non-simple points from list
            di3 = ~p_is_simple(nhood);
%             nhood(di3,:)=[];
%             cands(di3)=[];
            x(di3)=[];
            y(di3)=[];
            z(di3)=[];
            
            
            % if any candidates left: divide into 8 independent subvolumes
            if (~isempty(x))
                x1 = logical(mod(x,2));
                x2 = ~x1;
                y1 = logical(mod(y,2));
                y2 = ~y1;
                z1 = logical(mod(z,2));
                z2 = ~z1;
                ilst(1).l = x1 & y1 & z1;
                ilst(2).l = x2 & y1 & z1;
                ilst(3).l = x1 & y2 & z1;
                ilst(4).l = x2 & y2 & z1;
                ilst(5).l = x1 & y1 & z2;
                ilst(6).l = x2 & y1 & z2;
                ilst(7).l = x1 & y2 & z2;
                ilst(8).l = x2 & y2 & z2;
                
%                 idx = [];
                
                % do parallel re-checking for all points in each subvolume
                for i = 1:8                    
                    if any(ilst(i).l)
                        idx = ilst(i).l;
                        li = sub2ind([width height depth],x(idx),y(idx),z(idx));
                        skel(li)=0; % remove points
                        nh = pk_get_nh(skel,li);
                        di_rc = ~p_is_simple(nh);
                        if any(di_rc) % if topology changed: revert
                            skel(li(di_rc)) = true;
                        else
                            noChange = false; % at least one voxel removed
                        end
                    end
                end
            end
        end
        
        if( noChange )
            unchangedBorders = unchangedBorders + 1;
        end
        
    end
end

% get rid of padded zeros
skel = skel(2:end-1,2:end-1,2:end-1);
end


