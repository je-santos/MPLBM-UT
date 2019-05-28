function [vox,n_idx,ep] = pk_follow_link(skel,node,k,j,idx,cans,c2n)

vox = [];
n_idx = [];
ep = 0;

% assign start node to first voxel
vox(1) = node(k).idx(j);

i=1;
isdone = false;
while(~isdone) % while no node reached
    i=i+1; % next voxel
    next_cand = c2n(idx);
        cand = cans(next_cand,2);
        if(cand==vox(i-1)) % switch direction
            cand = cans(next_cand,3);
        end;
        if(skel(cand)>1) % node found
            vox(i) = idx;
            vox(i+1) = cand; % first node
            n_idx = skel(cand)-1; % node #
            if(node(n_idx).ep)
                ep=1;
            end;
            isdone = 1;
        else % next voxel
            vox(i) = idx;
            idx = cand;
        end;
end;
