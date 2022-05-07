function[] = ...
            output_tensors_fromSim(folder_loc, sim_folder,sim_size,save_to)


factor = containers.Map;
factor('1')  = 1.619499463;				
factor('2')  = 0.790754568;
factor('5')  = 0.29330618;
factor('10') = 0.128986101;
factor('20') = 0.051011612;

        
pressure = split(sim_folder,'_');
pressure = pressure{end};

conv_factor = factor(pressure);
        
file_names = dir( [folder_loc '/' sim_folder '/' sim_folder '_uz'] );


if isfile([save_to '/' sim_folder '_uz.mat'])
    fprintf('Simulation file already exists \n')
   return 
end


if isempty(file_names)
    fprintf('No sim was found \n')
    return
end

its = zeros(length(file_names),1) ;
for i = 1:length(file_names)
    its(i) = length(file_names(i).name) ;
end

[~,idx] = sort(its,'descend');
file_names = file_names(idx);
it = file_names(1).name;
it(1:7) = []; 
iteration=it;

if it == 500000
    return
end


geom_=0;

numPx=4;
numPy=6;
numPz=6;

lx=sim_size;
ly=sim_size;
lz=sim_size+2;

%++++++++++++++++++END USER INPUTS++++++++++++++%

subx = floor(lx/numPx);
suby = floor(ly/numPy);
subz = floor(lz/numPz);
xstart(1)=1;
xend(numPx)=lx;
for a=2:numPx
    xstart(a)=subx*(a-1)+1;
    xend(a-1)=subx*(a-1);
end
ystart(1)=1;
yend(numPy)=ly;
for a=2:numPy
    ystart(a)=suby*(a-1)+1;
    yend(a-1)=suby*(a-1);
end
zstart(1)=1;
zend(numPz)=lz;
for a=2:numPz
    zstart(a)=subz*(a-1)+1; %temp fix, zstart(a)=subz*(a-1)+zstart(1)+1;
    zend(a-1)=subz*(a-1);  %temp fix, zend(a-1)=subz*(a-1)+zstart(1);
end

imageVolUx=zeros(lx,ly,lz);
imageVolUy=zeros(lx,ly,lz);
imageVolUz=zeros(lx,ly,lz);
imageVolRho=zeros(lx,ly,lz);
for rank = 0:(numPx*numPy*numPz)-1
    filenameUx=[folder_loc '/' sim_folder '/' sim_folder '_ux/ux_' num2str(rank) '_' num2str(iteration)]; 
    fIdUx=fopen(filenameUx,'r');
    Ux=fscanf(fIdUx,'%f');
    fclose(fIdUx);
    filenameUy=[folder_loc '/' sim_folder '/' sim_folder '_uy/uy_' num2str(rank) '_' num2str(iteration)]; 
    fIdUy=fopen(filenameUy,'r');
    Uy=fscanf(fIdUy,'%f');
    fclose(fIdUy);
    filenameUz=[folder_loc '/' sim_folder '/' sim_folder '_uz/uz_' num2str(rank) '_' num2str(iteration)]; 
    fIdUz=fopen(filenameUz,'r');
    Uz=fscanf(fIdUz,'%f');
    fclose(fIdUz);
    filenameRho=[folder_loc '/' sim_folder '/' sim_folder '_rho/rho_' num2str(rank) '_' num2str(iteration)]; 
    fIdRho=fopen(filenameRho,'r');
    Rho=fscanf(fIdRho,'%f');
    fclose(fIdRho);
    count=0;
    %disp(rank);
    threadk=floor(rank/(numPx*numPy))+1;
    threadj=floor((rank-((threadk-1)*numPx*numPy))/numPx)+1;
    threadi=rank-((threadj-1)*numPx)-((threadk-1)*numPx*numPy)+1;

    sizeX=xend(threadi)-xstart(threadi)+1;
    sizeY=yend(threadj)-ystart(threadj)+1;
    sizeZ=zend(threadk)-zstart(threadk)+1;
    arrayUx=zeros(sizeX,sizeY,sizeZ);
    arrayUy=zeros(sizeX,sizeY,sizeZ);
    arrayUz=zeros(sizeX,sizeY,sizeZ);
    arrayRho=zeros(sizeX,sizeY,sizeZ);
    if(geom_==1)
        filenameGeom=[foldername '/' foldername '_GeometryOut/GeometryOut_' num2str(rank)]; 
        fIdGeom=fopen(filenameGeom,'r');
        Geom=fscanf(fIdGeom,'%f');
        fclose(fIdGeom);
        arrayGeom=zeros(sizeX,sizeY,sizeZ);
    end
    for k = 1:sizeZ
        for j = 1:sizeY
            for i = 1:sizeX
                count=count+1;
                arrayUx(i,j,k)=Ux(count);
                arrayUy(i,j,k)=Uy(count);
                arrayUz(i,j,k)=Uz(count);
                arrayRho(i,j,k)=Rho(count);
                if geom_==1
                    arrayGeom(i,j,k)=Geom(count);
                end
            end
        end
    end
    imageVolUx(xstart(threadi):xend(threadi),ystart(threadj):yend(threadj),zstart(threadk):zend(threadk))=arrayUx;
    imageVolUy(xstart(threadi):xend(threadi),ystart(threadj):yend(threadj),zstart(threadk):zend(threadk))=arrayUy;
    imageVolUz(xstart(threadi):xend(threadi),ystart(threadj):yend(threadj),zstart(threadk):zend(threadk))=arrayUz;
    imageVolRho(xstart(threadi):xend(threadi),ystart(threadj):yend(threadj),zstart(threadk):zend(threadk))=arrayRho;
    if geom_==1
        imageVolGeom(xstart(threadi):xend(threadi),ystart(threadj):yend(threadj),zstart(threadk):zend(threadk))=arrayGeom;
    end
    %imagesc(imageVolRho(:,:,zstart(threadk))),drawnow;
end

ux  = chop_bc( imageVolUx*conv_factor );
uy  = chop_bc( imageVolUy*conv_factor );
uz  = chop_bc( imageVolUz*conv_factor );
rho = chop_bc( imageVolRho );

save([save_to '/' sim_folder '_ux' ],'ux' );
save([save_to '/' sim_folder '_uy' ],'uy' );
save([save_to '/' sim_folder '_uz' ],'uz' );
save([save_to '/' sim_folder '_rho'],'rho');

mfp_folder = erase(folder_loc, '../domains');
mfp = chop_bc( unpackStruct( load([ '../mfp' mfp_folder '/' sim_folder '_MPa']) ) );

mfp(mfp==Inf) = 0;

save([save_to '/' sim_folder], 'mfp', 'ux', 'uy', 'uz', 'rho');

if strcmp(pressure, '20')
    
    conv_factor = ( (1e-6)/(0.5e-9) )^2;
    ux  = ux*conv_factor;
    uy  = uy*conv_factor;
    uz  = uz*conv_factor;
    
    lbm_name = erase(sim_folder, ['_' pressure]) ;
    save([save_to '/' lbm_name], 'ux', 'uy', 'uz', 'rho');
  
end
    
    
    

% 
% %calculate darcy uZ
% darcyUz=sum(sum(sum(imageVolUz)))/(lx*ly*lz);
% meanUz=sum(sum(sum(imageVolUz)))/sum(sum(sum(imageVolRho>0)));
% 
% for k = 1:lz
%     sliceMeanUz(k,1)=sum(sum(imageVolUz(:,:,k)))/sum(sum(imageVolRho(:,:,k)>0));
%     sliceMeanRho(k,1)=sum(sum(imageVolRho(:,:,k)))/sum(sum(imageVolRho(:,:,k)>0));
% end
% 
% for k = 1:lz
%     poreSum(k)=0;
%     for j = 1:ly
%         for i = 1:lx
%             if imageVolRho(i,j,k)>0 %&& imageVolUz(i,j,k)>0
%                 poreSum(k)=poreSum(k)+1;
%             end
%         end
%     end
%     porosityZ(k,1)=poreSum(k)/(lx*ly);
% end
% porosity=sum(poreSum)/(lx*ly*lz); %%????the porosity is 1.00000 ????
% fprintf('mean darcy: %f\n mean velocity: %f\n porosity: %f\n', darcyUz, meanUz, porosity);
% sliceUz=reshape(imageVolUz(:,:,1),lx,ly);
% maxUz=max(max(sliceUz));
% 
% % %-------------------------------------------------
% %Ux_P=imageVolUx*usurf_PL/(1e-8/accelz_L);
% %Uy_P=imageVolUy*usurf_PL/(1e-8/accelz_L);
% %Uz_P=imageVolUz*usurf_PL/(1e-8/accelz_L);
% % saveVTK_Vel(Ux_P,Uy_P,Uz_P,'Vel_P_20MPa');
% %------------------------------------------
% %save (['Uz_P_',num2str(int8(p0_P/1000000)),'MPa'], 'Uz_P');
% %save (['Vel_P_',num2str(int8(p0_P/1000000)),'MPa'], 'Ux_P','Uy_P','Uz_P');
end

function im = chop_bc(im)
    im(:,:,1)   = [];
    im(:,:,end) = [];
end