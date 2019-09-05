%% Read Data


f.f1_vti_struct=xml2struct('blob.vti'); %read output file
% %f.f1_vti_struct=xml2struct('3_tubes.vti'); %read output file
%f.f1_vti_struct=xml2struct('3tubes_pushback.vti'); %read output file
% %f.f2_vti_struct=xml2struct('rho2.vti'); %read output file
f.f1_vti_str=base64decode(f.f1_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
% %f.f2_vti_str=base64decode(f.f2_vti_struct.VTKFile.ImageData.Piece.PointData.DataArray.Text);
f.f1_vti_no=typecast([0 0 f.f1_vti_str],'double');
% %f.f2_vti_no=typecast([0 0 f.f2_vti_str],'double');
vti_size=str2num(f.f1_vti_struct.VTKFile.ImageData.Attributes.WholeExtent); %#ok<ST2NM>
vti_x=vti_size(1)+vti_size(2)+1;
vti_y=vti_size(3)+vti_size(4)+1;
vti_z=vti_size(5)+vti_size(6)+1;
rhof1=reshape(f.f1_vti_no(2:end),[vti_x vti_y vti_z]);
% %rhof2=reshape(f.f2_vti_no(2:end),[vti_x vti_y vti_z]);
clear f


tp_x=[];
tp_y=[];
tp_z=[];

wall_all=[];
angle_all=[];
angle_w=[];
angle_nw=[];
angle_n=[];

no_points=2;


%% Operations
for sl=100:110
    slice=squeeze(rhof1(:,:,sl)); %slice the 3d array
    image(slice*100);hold on;
    
    slice1=slice;
    slice1(slice==0 | slice==0.4 | slice==-0.4)=nan; %'remove' solid surface
    
    slice2=slice;
    slice2(slice==0 | slice==0.4 | slice==-0.4)=-100; %all the solid pixels n,nw
    
    [tt1]=contourc(slice1,[1 1]); %contour the meniscae
    if isempty(tt1)==0
        
        [x_f,y_f,z_f]=C2xyz(tt1); %arrange the results of the contour matrix
        clear tt1
        
        [tt2]=contourc(slice2,[-99 -99]); %contour the solid
        [x_s,y_s,z_s]=C2xyz(tt2);
        clear tt2
        
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
        
        scatter(x_sol,y_sol,'filled'); hold on;
        
        for mm=1:numel(x_f) %number of meniscae
            xftmp=x_f{mm}; %work with individual meniscous
            yftmp=y_f{mm};
            scatter(xftmp,yftmp,'r','filled'); hold on
            
            first_points_x=[xftmp(1),xftmp(end)];
            first_points_y=[yftmp(1),yftmp(end)];
            
            
            
            for i=1:numel(x_sol)
                p_dist(i)=sqrt((first_points_x(1)-x_sol(i))^2+(first_points_y(1)-y_sol(i))^2);
            end
            
            ind=find(p_dist==min(p_dist));
            ind=ind(1);
            triplepoint_x=x_sol(ind);
            triplepoint_y=y_sol(ind);
            
            tp_x(end+1)=triplepoint_x;
            tp_y(end+1)=triplepoint_y;
            tp_z(end+1)=sl;
            
            if (first_points_x(1)~= first_points_x(2) | first_points_y(1)~= first_points_y(2))==1
                
                
                scatter(triplepoint_x,triplepoint_y,100,'filled','k');hold on;
                line([xftmp(1+no_points),xftmp(1)],[yftmp(1+no_points),yftmp(1)]); hold on;
                line([x_sol(ind+1),x_sol(ind-1)],[y_sol(ind+1),y_sol(ind-1)]);
                
                
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
                
                uu=[(yftmp(1+no_points)-yftmp(1)),(xftmp(1+no_points)-xftmp(1)) ];
                angle1(mm)=acosd(dot(uu,vv)/(norm(uu)*norm(vv)));
                angle_all(end+1)=angle1(mm);
                wall_all(end+1)=slice(round(triplepoint_y),round(triplepoint_x));
                
                
                text(triplepoint_x+2,triplepoint_y+5,num2str(round(angle1(mm))),'FontSize',14,'BackgroundColor','w');
                
                
                %%%%%%%%second angle in meniscous
                
                for i=1:numel(x_sol)
                    p_dist(i)=sqrt((first_points_x(2)-x_sol(i))^2+(first_points_y(2)-y_sol(i))^2);
                end
                ind=find(p_dist==min(p_dist));
                ind=ind(1);
                triplepoint_x=x_sol(ind);
                triplepoint_y=y_sol(ind);
                
                tp_x(end+1)=triplepoint_x;
                tp_y(end+1)=triplepoint_y;
                tp_z(end+1)=sl;
                
                scatter(triplepoint_x,triplepoint_y,100,'filled','k');hold on;
                line([xftmp(end-no_points),xftmp(end)],[yftmp(end-no_points),yftmp(end)]); hold on;
                line([x_sol(ind+1),x_sol(ind-1)],[y_sol(ind+1),y_sol(ind-1)]);
                
                
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
                
                uu=[(yftmp(end-no_points)-yftmp(end)),(xftmp(end-no_points)-xftmp(end)) ];
                angle2(mm)=acosd(dot(uu,vv)/(norm(uu)*norm(vv)));
                angle_all(end+1)=angle2(mm);
                wall_all(end+1)=slice(round(triplepoint_y),round(triplepoint_x));
                text(triplepoint_x-3,triplepoint_y+5,num2str(round(angle2(mm))),'FontSize',14,'BackgroundColor','w');
                clear vv
                
            else
                angle1(mm)=180;
                angle2(mm)=180;
                text(triplepoint_x-3,triplepoint_y+5,num2str(round(angle2(mm))),'FontSize',14,'BackgroundColor','w');
            end
            
        end
        
        %pause(.5)
    end
end

figure;scatter(wall_all,angle_all);
xlabel('<---------------- Water Wet');
title('Contact Angle Sensitivity Analysis');
ylabel('Angle [degrees]');