function perc_length = Test_Skel2Graph3D(skel,testvol)
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




counter=1;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
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
    NODE(i,3)=node2(i).comx;
    NODE(i,4)=node2(i).comx;
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
        
        for k=1:count           
            [x3,y3,z3]=ind2sub([w,l,h],link2(j).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(j).point(k+1));
            if k==1
              P1 = [x3,y3,z3]; 
            end
            if k==count
              P2 = [x2,y2,z2];  
            end
            D = sqrt((x2-x3)^2+(y2-y3)^2+(z2-z3)^2)+D;
        end
       link2(j).dist=D;
       link2(j).point1=P1;
       link2(j).point2=P2;
       
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

LINK(j,1)=j;
LINK(j,2)=link2(j).n1;
LINK(j,3)=link2(j).n2;
LINK(j,4)=D;
       
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
        

        
        
%         path1(p,2)=x3;
%         path1(p,3)=y3;
%         path1(p,4)=z3;
%         
%         path1(p,5)=x2;
%         path1(p,6)=y2;
%         path1(p,7)=z2;
        
%         [x3,y3,z3]=ind2sub([w,l,h],link2(path(p)).point(k));
%         [x2,y2,z2]=ind2sub([w,l,h],link2(path(p+1)).point(k));

 %       line([x3 x2],[y3 y2],[z3 z2],'Color',col,'LineWidth',2);


end



% Inputs:
%     NODES should be an Nx3 or Nx4 matrix with the format [ID X Y] or [ID X Y Z]
%       where ID is an integer, and X, Y, Z are cartesian position coordinates)
%     SEGMENTS should be an Mx3 matrix with the format [ID N1 N2 D]
%       where ID is an integer, and N1, N2 correspond to node IDs from NODES list
%       such that there is an [undirected] edge/segment between node N1 and
%       node N2, and D is calculated tortuous distance between nodes
%     SID should be an integer in the node ID list corresponding with the starting node
%     FID (optional) should be an integer in the node ID list corresponding with the finish

 
% numb=numel(pores_in);
% for tt=1:numb
%     P=pores_in(tt);
%     PP(tt)=Pores(P).rInscribed;
% end
% maxPP =max(PP);
% for tt=1:numb
% if PP(tt) == maxPP
% PoreI=pores_in(tt);
% end
% end
% % PoreF=PoreI;
% % PoreI=Finalpore;



% unvisited=zeros([1 npores]);
% stoploop=2;
% counter=0;
% D_PoreI=0;
% Dist= Inf([1,npores]);
% Dist(PoreI)=0;
% while PoreI~=Finalpore 
% counter=counter+1;
% 
% P_neighbors=Pores(PoreI).adjPores;
% T_neighbors=Pores(PoreI).adjThroats;
% N1=numel(P_neighbors);
% Ci = Pores(PoreI).center;
% 
% Dij=0;
% 
% for tt=1:N1
%     P_N=P_neighbors(tt);
%     T_N=T_neighbors(tt);
% %    if indicatorPores(P_N)== flag_nonwet
%      if P_N~=PoreI
%          if unvisited(P_N)== 0
%         Cj =  Pores(P_N).center;
%          else
%         Cj = Ci;
%          end
%         Dij(tt) = ((Cj(1)-Ci(1))^2+(Cj(2)-Ci(2))^2+(Cj(3)-Ci(3))^2)^0.5+D_PoreI;
%     
%        if Dij(tt)<Dist(P_N)
%         Dist(P_N)=Dij(tt);
%         PreviousP(P_N)=PoreI;
% %       PreviousT(T_N)=PoreI;        
%        end
%     else
%         Dij(tt)=D_PoreI;
%     end
% end
% unvisited(PoreI)=1;
% %unvisited(PoreI)=1;
% 
% %MiniD=min(Dij(Dij~=D_PoreI)); % find minimum D from source for all unvisited =0 pores and designate next PoreI
% MiniD=min(Dist(unvisited==0));
% 
% 
% for tt=1:npores
%     if MiniD==Dist(tt)&&unvisited(tt)==0
% %       D=D+Dij(tt);
% %       PoreOld1=PoreOld;
% %       PoreOld=PoreI;
%         PoreI=tt;
%         D_PoreI=MiniD;
%         break
%     end
% end
%   
% %     check = indicatorPores(PoreI)*xlo_side(PoreI); 
% %  if check == 16
% %      stoploop =0;    
% %  end
% end
% 


perc_length=dist;

