function [NODE,LINK,perc_length] = Test_Skel2Graph3D(skel,testvol)
% clear all;
% close all;

% load example binary skeleton image
% load skel

% load example binary skeleton image
% load skel_test


figure();
col=[.7 .7 .8];
hiso = patch(isosurface(testvol,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(testvol,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(testvol,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
% w=size(skel,1);
% l=size(skel,2);
% h=size(skel,3);
% [x,y,z]=ind2sub([w,l,h],find(skel(:)));
% plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');
% set(gcf,'Color','white');
view(140,80)


w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node,link] = Skel2Graph3D(skel,0);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,0);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)
    
    wl = wl_new;
    
    skel2 = Graph2Skel3D(node2,link2,w,l,h);
    [A2,node2,link2] = Skel2Graph3D(skel2,0);
    
    wl_new = sum(cellfun('length',{node2.links}));
    
end

distance = bwdist(~testvol);
% distance = skel.*dist_trans;


counter=1;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    x_1=round(x1);
    y_1=round(y1);
    z_1=round(z1);
    
    poreR(i)= distance(x_1,y_1,z_1);
    
    X1(i)=x1;
    
    [Xmax,Imax] = max(X1);
    [xmin,imin] = min(X1);
    
    if X1(i) < 12
        X_minoption(i) = x1;
        Y_minoption(i) = y1;
        Z_minoption(i) = z1;
    else
        X_minoption(i) = inf;
        Y_minoption(i) = inf;
        Z_minoption(i) = inf;
        
    end
    
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'c';
    end
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='k'; % links are red
        end
        if(node2(i).ep==1)
            col='k';
        end
        
        
        % draw edges as lines using voxel positions
        count=length(link2(node2(i).links(j)).point)-1 ;
        for k=1:count
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
            
            
        end
        
    end
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end

for i=1:length(node2)
    Ymax = node2(Imax).comy;
    Zmax = node2(Imax).comz;
    if Z_minoption(i) > 0 && Z_minoption(i) < inf
        Distance(i)= sqrt((Y_minoption(i)-Ymax)^2+(Z_minoption(i)-Zmax)^2);
    else
        Distance(i)= inf;
    end
    NODE(i,1)=i;
    NODE(i,2)=node2(i).comx;
    NODE(i,3)=node2(i).comy;
    NODE(i,4)=node2(i).comz;
    NODE(i,5)= poreR(i);
end
[Dmin,Imin] = min(Distance);
Xmin = node2(Imin).comx;
Ymin = node2(Imin).comy;
Zmin = node2(Imin).comz;





for j=1:numel(link2)
    % if Pointer1(count,1)== node2(link2(count).n1).comx && Pointer1(count,2)== node2(link2(count).n1).comy && Pointer1(count,3)== node2(link2(count).n1).comz
    D=0;
    % draw edges as lines using voxel positions
    % k=node2(i).links(j);
    count=length(link2(j).point)-1 ;
    throat=0;
    for k=1:count
        [x3,y3,z3]=ind2sub([w,l,h],link2(j).point(k));
        [x2,y2,z2]=ind2sub([w,l,h],link2(j).point(k+1));
        throat(k)= distance(x3,y3,z3);
        if k==1
            P1 = [x3,y3,z3];
        end
        if k==count
            P2 = [x2,y2,z2];
            throat(count+1)= distance(x2,y2,z2);
        end
        D = sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)+D;
    end
    link2(j).dist=D;
    link2(j).point1=P1;
    link2(j).point2=P2;
    link2(j).throatR=min(throat);
    
    %        Pointer1X = node2(link2(j).n1).comx;
    %        Pointer1Y = node2(link2(j).n1).comy;
    %        Pointer1Z = node2(link2(j).n1).comz;
    %        Pointer2X = node2(link2(j).n2).comx;
    %        Pointer2Y = node2(link2(j).n2).comy;
    %        Pointer2Z = node2(link2(j).n2).comz;
    %
    %        Pointer1=[Pointer1X Pointer1Y Pointer1Z];
    %        Pointer2=[Pointer2X Pointer2Y Pointer2Z];
    %
    %        link2(j).pointer1=Pointer1;
    %        link2(j).pointer2=Pointer2;
    
    %throat(throat==0)=inf;
    
    LINK(j,1)=j;
    LINK(j,2)=link2(j).n1;
    LINK(j,3)=link2(j).n2;
    LINK(j,4)=D;
    LINK(j,5)= min(throat);
    LINK(j,6)= mean(throat);
end

% Distance=unique(Distance,'stable');
% Pointer1=unique(Pointer1,'stable');
% Pointer2=unique(Pointer2,'stable');
% link2.dist=Distance;
% link2.point1=Pointer1;
% link2.point2=Pointer2;

axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);
pliiist=0;

% % Find single shortest path  by Dijkstras algorithm (works)
start_id = Imax;
finish_id = Imin;
perc_length=0;

DJ=1;

if DJ==1
    
    disp('Finding shortest path  by Dijkstras algorithm');
    [dist,path] = dijkstra(NODE,LINK,start_id,finish_id);
    path=path';
    
    for p=1:(numel(path)-1)
        
        col='r'; % links are red
        % draw edges as lines using voxel positions
        
        %     count=length(link2(path(p)).point)-1 ;
        %     for k=1:count
        %         [x3,y3,z3]=ind2sub([w,l,h],link2(path(p)).point(k));
        %         [x2,y2,z2]=ind2sub([w,l,h],link2(path(p)).point(k+1));
        %         line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        %     end
        %         [x3,y3,z3]=ind2sub([w,l,h],link2(path(p)).point(k));
        %         [x2,y2,z2]=ind2sub([w,l,h],link2(path(p+1)).point(k));
        %         line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        
        x3=node2(path(p)).comx;
        y3=node2(path(p)).comy;
        z3=node2(path(p)).comz;
        
        x2=node2(path(p+1)).comx;
        y2=node2(path(p+1)).comy;
        z2=node2(path(p+1)).comz;
        
        P1=[x3 y3 z3];
        P2=[x2 y2 z2];
        
        num=numel(link2);
        
        %      line([x3 x2],[y3 y2],[z3 z2],'Color',col,'LineWidth',2);
        
        N1=path(p);
        N2=path(p+1);
        
        
        
        for count=1:num
            NN1= link2(count).n1;
            NN2= link2(count).n2;
            if(NN1 == N1 && NN2 == N2 || NN1 == N2 && NN2 == N1 )
                counter=length(link2(count).point)-1;
                for k=1:counter
                    [x3,y3,z3]=ind2sub([w,l,h],link2(count).point(k));
                    [x2,y2,z2]=ind2sub([w,l,h],link2(count).point(k+1));
                    line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
                end
                break
            end
        end
    end
    perc_length=dist;
end

