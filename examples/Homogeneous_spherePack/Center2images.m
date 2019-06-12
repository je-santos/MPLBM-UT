%% convert sphere pack with center coordinates and radii to image slices
%Abhishek Bihani, March 2019

% Input Data can be found here: https://www.digitalrocksportal.org/projects/204

tic
clear
data = load('C:\Users\Abhishek\Desktop\16th dec 18 work\gravity_final\excels\10to1_50_1.csv');

x = data(:,1);  % x-coordinate
y = data(:,2);  % y-coordinate
z = data(:,3);  % z-coordinate
r = data(:,4); % radius

a = min(r);
b = max(r);
a1 = nnz(r==a);
b1 = nnz(r==b);
Vl= b1*4/3*pi*b^3;
Vs= a1*4/3*pi*a^3;
VL= Vl*100/(Vs+Vl);
VL = round(VL,2);  % Volume percent of large grains
ratio = b/a;       % Ratio of radius of large to small grains
ratio=round(ratio);

filename = sprintf('sliced_grainpack_%gto1_VL%g_',ratio,VL);% create output name

numb=length(x);

for ii=1:numb
    x(ii)=x(ii)*2*10^4;
    y(ii)=y(ii)*2*10^4;
    z(ii)=z(ii)*2*10^4;
    r(ii)=r(ii)*10^4;
end

res=50;  % resolution

% Finding max and min of grainpack
[Xmax,I]= max(x);
Xmax=Xmax+r(I);
[Xmin,I]= min(x);
Xmin=Xmin-r(I);
[Ymax,I]= max(y);
Ymax=Ymax+r(I);
[Ymin,I]= min(y);
Ymin=Ymin-r(I);
[Zmax,I]= max(z);
Zmax=Zmax+r(I);
[Zmin,I]= min(z);
Zmin=Zmin-r(I);
% Xmin=0;
% Ymin=0;
% Zmin=0;
X_len= round((Xmax-Xmin)*res);
Y_len= round((Ymax-Ymin)*res);
Z_len= round((Zmax-Zmin)*res);

D_len=[X_len,Y_len,Z_len];
D_len=max(D_len);

%bw=zeros(X_len,Y_len,Z_len);

lims = [1 D_len];
[x_mesh,y_mesh,z_mesh] = meshgrid(lims(1):lims(2));

%[x_mesh,y_mesh,z_mesh]=meshgrid(X_len,Y_len,Z_len);

%Populating grains using meshgrid

for i=1:numb
    X=x(i)*res;
    Y=y(i)*res;
    Z=z(i)*res;
    R=r(i)*res;
    
    bw1=sqrt((x_mesh-X).^2+(y_mesh-Y).^2+(z_mesh-Z).^2)<=R;
    if i==1
        bw=bw1;
    else
        bw=bw1|bw;
    end
end

simSpace=bw;

write_images=0;

if write_images ==1
    % Writing image slices
    mkdir(filename);
    sizeK=size(simSpace,3);
    for k = 1:sizeK
        k1=sprintf('%g.bmp',k);
        imwrite(simSpace(:,:,k),[filename '/' filename num2str(k) '.bmp'],'bmp');
    end
end

toc

filename1=[filename '.mat'];
save(filename1,'simSpace');

addpath ('../../pre-processing') %pre-processing libraries
[status]=create_LBM_MAT2DAT(simSpace,filename); % Creating Palabos input file

