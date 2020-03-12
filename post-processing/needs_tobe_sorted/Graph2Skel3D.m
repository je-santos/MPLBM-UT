function skel = Graph2Skel3D(node,link,w,l,h)

% create binary image
skel = false(w,l,h);

% for all nodes
for i=1:length(node)
    if(~isempty(node(i).links)) % if node has links
        skel(node(i).idx)=true; % node voxels
        a = [link(node(i).links(node(i).links>0)).point];
        if(~isempty(a))
            skel(a)=1; % edge voxels
        end;
    end;
end;

