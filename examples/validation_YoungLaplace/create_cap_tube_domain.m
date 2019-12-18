%This script creates a number of parallel capillary tubes of specified
%radius to test the capillary pressure needed to flow through each tube

radius=[31,21,16,11,6,5,4,3];

no_circles=1;
spacing_between=5; %spacing between adjacent tubes
mesh=1; % if 1x1 mesh wanted near outlet

radius_c = max(radius);
dom.size_x=max(radius_c*2)+spacing_between*2;
dom.size_y=sum(radius_c*2)+spacing_between*(numel(radius_c)+1);
dom.size_z=50; %length

B=zeros(dom.size_x,dom.size_y,dom.size_z); %3D matrix
C=zeros(dom.size_x,dom.size_y);  %fracture walls
D=ones(dom.size_x,dom.size_y)*2;  %fracture walls
DDD=ones(dom.size_x,dom.size_y)*2;  %fracture walls
J=ones(dom.size_x,dom.size_y)*2;

dom.radius=radius_c;
dom.value=[1,1,1];
dom.xcenters=[dom.size_x/2];

last=0;
for i=1:numel(radius_c)
    dom.ycenters(i)=[spacing_between+radius_c(i)]+last;
    last=dom.ycenters(i)+radius_c(i);
end

for i=1:numel(dom.ycenters)
    
    xc = int16(dom.xcenters);
    yc = int16(dom.ycenters(i));
    
    x = int16(0);
    y = int16(dom.radius(i));
    d = int16(1 - dom.radius);
    
    C(xc, yc+y) = dom.value(i);
    C(xc, yc-y) = dom.value(i);
    C(xc+y, yc) = dom.value(i);
    C(xc-y, yc) = dom.value(i);
    
    while ( x < y - 1 )
        x = x + 1;
        if ( d < 0 )
            d = d + x + x + 1;
        else
            y = y - 1;
            a = x - y + 1;
            d = d + a + a;
        end
        C( x+xc,  y+yc) = dom.value(i);
        C( y+xc,  x+yc) = dom.value(i);
        C( y+xc, -x+yc) = dom.value(i);
        C( x+xc, -y+yc) = dom.value(i);
        C(-x+xc, -y+yc) = dom.value(i);
        C(-y+xc, -x+yc) = dom.value(i);
        C(-y+xc,  x+yc) = dom.value(i);
        C(-x+xc,  y+yc) = dom.value(i);
    end
    
    
    for ii = xc-int16(dom.radius(i)):xc+(int16(dom.radius(i)))
        for jj = yc-int16(dom.radius(i)):yc+(int16(dom.radius(i)))
            tempR = sqrt((double(ii) - double(xc)).^2 + (double(jj) - double(yc)).^2);
            if(tempR <= double(int16(dom.radius(i))))
                D(ii,jj)=0;
            end
        end
    end
    
    D(C(:,:)==dom.value(i))=dom.value(i);
end

ttubes(:,:,1)=D;
first_tube=D;
D=D*0;
D=D+2;
C=C*0;


%%%%%%%%%%%%%%%%%%%%% Second loop
for t=1:(numel(radius)-1)
    dom.radius=radius(t+1);
    
    for i=1:numel(dom.ycenters)
        
        xc = int16(dom.xcenters);
        yc = int16(dom.ycenters(i));
        
        x = int16(0);
        y = int16(dom.radius(i));
        d = int16(1 - dom.radius);
        
        C(xc, yc+y) = dom.value(i);
        C(xc, yc-y) = dom.value(i);
        C(xc+y, yc) = dom.value(i);
        C(xc-y, yc) = dom.value(i);
        
        while ( x < y - 1 )
            x = x + 1;
            if ( d < 0 )
                d = d + x + x + 1;
            else
                y = y - 1;
                a = x - y + 1;
                d = d + a + a;
            end
            C( x+xc,  y+yc) = dom.value(i);
            C( y+xc,  x+yc) = dom.value(i);
            C( y+xc, -x+yc) = dom.value(i);
            C( x+xc, -y+yc) = dom.value(i);
            C(-x+xc, -y+yc) = dom.value(i);
            C(-y+xc, -x+yc) = dom.value(i);
            C(-y+xc,  x+yc) = dom.value(i);
            C(-x+xc,  y+yc) = dom.value(i);
        end
        
        
        for ii = xc-int16(dom.radius(i)):xc+(int16(dom.radius(i)))
            for jj = yc-int16(dom.radius(i)):yc+(int16(dom.radius(i)))
                tempR = sqrt((double(ii) - double(xc)).^2 + (double(jj) - double(yc)).^2);
                if(tempR <= double(int16(dom.radius(i))))
                    D(ii,jj)=0;
                end
            end
        end
        
        D(C(:,:)==dom.value(i))=dom.value(i);
    end
    
    ttubes(:,:,t+1)=D;
    layers_tmp=ttubes(:,:,t);
    layers_tmp(layers_tmp==0)=1;
    layers_tmp(ttubes(:,:,t+1)==0)=0;
    
    layers(:,:,t)=layers_tmp;
    D=D*0;
    D=D+2;
    C=C*0;
        
end


domain=[];
%% print
final=[];
for i=1:numel(radius)
    tmp=[];
    tmp1=[];
    tmp=ttubes(:,:,i);
    %tmp1=tmp([36-radius(i)-4:36+radius(i)+4],:);
    tmp1=tmp([dom.size_x/2-radius(i)-4:dom.size_x/2+radius(i)+4],:);
    image(tmp1,'CDataMapping','scaled');
    axis equal
    [aa,bb]=size(tmp1);
    final=[final;tmp1];
    close
    %if i==1
    %final=tmp1;
    %else
    %   final([end+1:bb],:)=tmp1;
    %end
end
[r,c]=size(final);
c=c+20;
final=[2*ones(r,10) final 2*ones(r,10)];
final=[2*ones(10,c) ;  final ; 2*ones(10,c)];

[r,c]=size(final);

%NewFileName = sprintf('tubes_%g_%g_%g',r,c,dom.size_z);% create output directory
NewFileName = sprintf('tubes_%g_%g_%g',dom.size_z,c,r);% create output directory
NewFileName = [NewFileName '.dat'];

fid_domain1=fopen(['input/' NewFileName],'w');


for j=1:1
    fprintf(fid_domain1,'%i\n',final*0);
    domain(:,:,end)=final*0;
end

for j=1:1
    fprintf(fid_domain1,'%i\n',final*0);
    domain(:,:,end+1)=final*0;
end

for j=1:1
    tmp = final;
    tmp(tmp==2)=1;
    fprintf(fid_domain1,'%i\n',tmp);
    domain(:,:,end+1)=tmp;
end


for j=1:dom.size_z-6
    fprintf(fid_domain1,'%i\n',final);
    domain(:,:,end+1)=final;
end

for j=1:1
    tmp = final;
    tmp(tmp==2)=1;
    fprintf(fid_domain1,'%i\n',tmp);
    domain(:,:,end+1)=final*0;
end

if mesh ==1
    for j=1:1
        tmp1=toeplitz(mod(1:r,2),mod(1:c,2));
        tmp1(tmp1==1)=4;
        fprintf(fid_domain1,'%i\n',tmp1);
        domain(:,:,end+1)=tmp1;
    end
end

for j=1:1
    fprintf(fid_domain1,'%i\n',final*0);
    domain(:,:,end+1)=final*0;
end



fclose(fid_domain1);
%isosurface(B)

figure()

imagesc(final);
axis equal;

disp('2. Surface Created')