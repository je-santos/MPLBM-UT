%% Read Data
f.f1_vti_struct=xml2struct('rho1.vti'); %read output file
%f.f2_vti_struct=xml2struct('rho2.vti'); %read output file
f.f1_vti_str=base64decode(f.f1_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
%f.f2_vti_str=base64decode(f.f2_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
f.f1_vti_no=typecast([0 0 f.f1_vti_str],'double');
%f.f2_vti_no=typecast([0 0 f.f2_vti_str],'double');
vti_size=str2num(f.f1_vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
vti_x=vti_size(1)+vti_size(2)+1;
vti_y=vti_size(3)+vti_size(4)+1;
vti_z=vti_size(5)+vti_size(6)+1;
rhof1=reshape(f.f1_vti_no(2:end),[vti_x vti_y vti_z]);
%rhof2=reshape(f.f2_vti_no(2:end),[vti_x vti_y vti_z]);
clear f


%% Operations

%------------------------flatten down the surface
%------------------------linearly interpolate fluid and surface

slice=rhof1(:,:,105); %slice the 3d array

slice1=slice;
slice1(slice<=0)=nan; %solid surface

slice2=slice;
slice2(slice==0 | slice==0.4 | slice==-0.4)=-100; %all the solid pixels n,nw



[tt1]=contourc(slice1,[1 1]); %contour the meniscae
[x_f,y_f,z_f]=C2xyz(tt1); %arrange the results of the contour matrix



for mm=1:numel(x_f) %number of meniscae
    xftmp=x_f{mm};
    yftmp=y_f{mm};
    scatter(xftmp,yftmp,'r'); hold on
    
    
    [tt2]=contourc(slice2,[-99 -99]); %contour the solid
    [x_s,y_s,z_s]=C2xyz(tt2);
    
    
    
    xstmp=x_s{mm};
    ystmp=y_s{mm};
    
    scatter(xstmp,ystmp,'k');
    
    %put all the solids in an array
    k=1;
    for i=1:numel(x_s)
        clear x_tmp;
        clear y_tmp;
        x_tmp=double(x_s{i});
        y_tmp=double(y_s{i});
        for j=1:numel(x_tmp)
            x_sol(k)=x_tmp(j);
            y_sol(k)=y_tmp(j);
            k=k+1;
        end
    end
    %image(slice2*100);hold on;
    %scatter(x_sol,y_sol);
    
    first_points_x=[xftmp(1),xftmp(end)];
    first_points_y=[yftmp(1),yftmp(end)];
    
    for i=1:numel(x_sol)
        p_dist(i)=sqrt((first_points_x(1)-x_sol(i))^2+(first_points_y(1)-y_sol(i))^2);
    end
    
    ind=find(p_dist==min(p_dist));
    triplepoint_x=x_sol(ind);
    triplepoint_y=y_sol(ind);
    scatter(triplepoint_x,triplepoint_y,'filled');hold on;
    line([xftmp(3),xftmp(1)],[yftmp(3),yftmp(1)]); hold on;
    line([x_sol(ind+1),x_sol(ind-1)],[y_sol(ind+1),y_sol(ind-1)]);
    
    if mm==1
        mm
    end
    
    if slice2(round(triplepoint_y),round(triplepoint_x+1))>0
        if slice2(round(triplepoint_y+1),round(triplepoint_x+1))>slice2(round(triplepoint_y-1),round(triplepoint_x+1))
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind-1)-y_sol(ind+1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            end
        else
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind-1)+y_sol(ind+1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            end
            
        end
    elseif slice2(round(triplepoint_y),round(triplepoint_x-1))>0 %up and down
        if slice2(round(triplepoint_y+1),round(triplepoint_x-1))>slice2(round(triplepoint_y-1),round(triplepoint_x-1))
            if y_sol(ind+1)>y_sol(ind-1) %up and down
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind-1)-y_sol(ind+1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            end
        else
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind-1)+y_sol(ind+1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            end
            
        end
        
    elseif slice2(round(triplepoint_y+1),round(triplepoint_x))>0
        
        if slice2(round(triplepoint_y+1),round(triplepoint_x+1))>slice2(round(triplepoint_y+1),round(triplepoint_x-1))
            if x_sol(ind+1)>x_sol(ind-1) %left right
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind-1)-x_sol(ind+1))]; %solid normal
            end
        else
            if x_sol(ind+1)>x_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind-1)+x_sol(ind+1))]; %solid normal
            end
            
        end
        
    elseif slice2(round(triplepoint_y-1),round(triplepoint_x))>0
        if slice2(round(triplepoint_y-1),round(triplepoint_x+1))>slice2(round(triplepoint_y-1),round(triplepoint_x-1))
            if x_sol(ind+1)>x_sol(ind-1) %left right
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind-1)-x_sol(ind+1))]; %solid normal
            end
        else
            if x_sol(ind+1)>x_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind-1)+x_sol(ind+1))]; %solid normal
            end
            
        end
        
    else
        disp('Error')
        
    end
    
    
    u=[(xftmp(3)-xftmp(1)) , (yftmp(3)-yftmp(1))];
    uu=[(yftmp(3)-yftmp(1)),(xftmp(3)-xftmp(1)) ];
    v=[(x_sol(ind+1)-x_sol(ind-1)) , (y_sol(ind+1)-y_sol(ind-1))];
    angle=acosd(dot(u,v)/(norm(u)*norm(v)));
    angle1=acosd(dot(uu,vv)/(norm(uu)*norm(vv)))
    text(triplepoint_x+2,triplepoint_y+2,num2str(angle1));
    clear vv
    %second
    
    for i=1:numel(x_sol)
        p_dist(i)=sqrt((first_points_x(2)-x_sol(i))^2+(first_points_y(2)-y_sol(i))^2);
    end
    ind=find(p_dist==min(p_dist));
    triplepoint_x=x_sol(ind);
    triplepoint_y=y_sol(ind);
    scatter(triplepoint_x,triplepoint_y,'filled');hold on;
    line([xftmp(end-2),xftmp(end)],[yftmp(end-2),yftmp(end)]); hold on;
    line([x_sol(ind+1),x_sol(ind-1)],[y_sol(ind+1),y_sol(ind-1)]);
    
    
    
    % if slice2(round(triplepoint_y),round(triplepoint_x-1))>0
    %     if slice2(round(triplepoint_y+1),round(triplepoint_x-1))>slice2(round(triplepoint_y-1),round(triplepoint_x-1))
    %         if y_sol(ind+1)>y_sol(ind-1)
    %             vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
    %         else
    %             vv=[ (y_sol(ind-1)-y_sol(ind+1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
    %         end
    %     else
    %         if y_sol(ind+1)>y_sol(ind-1)
    %             vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
    %         else
    %             vv=[ (-y_sol(ind-1)+y_sol(ind+1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
    %         end
    %     end
    
    if slice2(round(triplepoint_y),round(triplepoint_x+1))>0
        if slice2(round(triplepoint_y+1),round(triplepoint_x+1))>slice2(round(triplepoint_y-1),round(triplepoint_x+1))
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind-1)-y_sol(ind+1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            end
        else
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind-1)+y_sol(ind+1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            end
            
        end
    elseif slice2(round(triplepoint_y),round(triplepoint_x-1))>0 %up and down
        if slice2(round(triplepoint_y+1),round(triplepoint_x-1))>slice2(round(triplepoint_y-1),round(triplepoint_x-1))
            if y_sol(ind+1)>y_sol(ind-1) %up and down
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind-1)-y_sol(ind+1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            end
        else
            if y_sol(ind+1)>y_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind-1)+y_sol(ind+1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            end
            
        end
        
    elseif slice2(round(triplepoint_y+1),round(triplepoint_x))>0
        
        if slice2(round(triplepoint_y+1),round(triplepoint_x+1))>slice2(round(triplepoint_y+1),round(triplepoint_x-1))
            if x_sol(ind+1)>x_sol(ind-1) %left right
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind-1)-x_sol(ind+1))]; %solid normal
            end
        else
            if x_sol(ind+1)>x_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind-1)+x_sol(ind+1))]; %solid normal
            end
            
        end
        
    elseif slice2(round(triplepoint_y-1),round(triplepoint_x))>0
        if slice2(round(triplepoint_y-1),round(triplepoint_x+1))>slice2(round(triplepoint_y-1),round(triplepoint_x-1))
            if x_sol(ind+1)>x_sol(ind-1) %left right
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind+1)-x_sol(ind-1))]; %solid normal
            else
                vv=[ (y_sol(ind+1)-y_sol(ind-1)),(x_sol(ind-1)-x_sol(ind+1))]; %solid normal
            end
        else
            if x_sol(ind+1)>x_sol(ind-1)
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind+1)+x_sol(ind-1))]; %solid normal
            else
                vv=[ (-y_sol(ind+1)+y_sol(ind-1)),(-x_sol(ind-1)+x_sol(ind+1))]; %solid normal
            end
            
        end
        
    else
        disp('Error')
        
    
end



u=[(xftmp(end-2)-xftmp(end)) , (yftmp(end-2)-yftmp(end))];
uu=[(yftmp(end-2)-yftmp(end)),(xftmp(end-2)-xftmp(end)) ];
v=[(x_sol(ind-1)-x_sol(ind+1)) , (y_sol(ind-1)-y_sol(ind+1))];
angle1=acosd(dot(uu,vv)/(norm(uu)*norm(vv)))
angle=acosd(dot(u,v)/(norm(u)*norm(v)));
text(triplepoint_x+2,triplepoint_y+2,num2str(angle1));
clear vv

end

xstmp=x_s{6};
ystmp=y_s{6};

scatter(xstmp,ystmp,'k');

