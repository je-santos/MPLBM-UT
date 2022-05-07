function [] = domains_4sim(sample_name, tiff_dir, mfp_loc, sim_geom_loc)

%sample_name = '10_01_256';
%tiff_dir = '../matlab_volumes_real_tiff';
%mfp_loc = '../mfp_real';
%sim_geom_loc = '../domains_real';


mfps      = [17.56, 8.7228, 3.4318, 1.6871, 0.84688];
pressures = [1,     2,      5,      10,     20];


periodicX_= 0;   %%'1' represents periodic boundary
periodicY_= 0;   %%
periodicZ_= 1;   %% 1' represents periodic boundary


addwallx_= 1;    %% 1 for porous media if the x is not periodic boundary
addwally_= 1;    %%

globalVisc_= 0;  %% set globalVisc_ to 1 to turn off local effective viscosity

sq=1.1; %dummy
a=1.00; %these are the alpha/beta values chosen in Landry et al. 2016 arctan function
beta=1.00;  %

% make sure that all are even and greater than 1
numPx=4;
numPy=6;
numPz=6*2;



for pres_count=1:length(pressures)
    
    pressure = pressures(pres_count)
    mfp = mfps(pres_count)
    
    
    
    im_tmp = imread([tiff_dir '/' sample_name '/' sample_name '_0001']);
    sim_size = size(im_tmp,1);
    
    lx=sim_size;         
    ly=sim_size;         
    sliceStart=1;
    sliceEnd=sim_size+2;   
    lz=sliceEnd-sliceStart+1;
    
    
    
    %%for simple tube, maxRad is diagonal of square, sqrt(2)*lx, why?
    maxRad=21;%ceil((2^0.5)*radTube);%for local mfp measurement, limits size of vol sent to svsqLoc3D
    %%for porous media with periodic structure, maxRad should be the diameter of the maximum pores
    
    if globalVisc_==1
        fileOutRoot=[sample_name '_' num2str(numPx) '-' num2str(numPy) '-' num2str(numPz) '_Global'];
    else
        fileOutRoot=[sample_name '_' num2str(numPx) '-' num2str(numPy) '-' num2str(numPz) '_' num2str(pressure) '_Psi_mfp'];% num2str(uint16(floor(mfp))) '-' num2str(uint16((mfp-floor(mfp))*1000))];
    end
    
    mkdir([sim_geom_loc '/' fileOutRoot]);
    
    
    cx  = [1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0];
    cy  = [0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1];
    cz  = [0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1];
    
    %%C1,C2 second-order slip boundary, Us=4*C1*Kn+8C2*Kn**2, Ref. Chai et al. 2010.Commun. Comput. Phys. 8(5), 1052–1073
    C1=1.11;
    C2=0.61;
    
    
    fluidNTotAll=0;     %%
    svTotAll=0;         %%
    sliceCount=0;
    
    if periodicY_==1
        svLocVol=zeros(lx+2,ly+2,lz); %%
        rLocVol=zeros(lx+2,ly+2,lz);
    else
        svLocVol=zeros(lx+2,ly+2,lz);
        rLocVol=zeros(lx+2,ly+2,lz);
    end
    outSlice=[];
    if globalVisc_==0
        %For LEV-LBM the sv and sq for every fluid node needs to be
        %determined following Landry et al. 2016, in the following each
        %slice in z is evaluated to build a volume containing the local sv
        %and sq values
        for slice = sliceStart:sliceEnd
            sliceCount=sliceCount+1;
            sliceLow=floor(slice-maxRad); %bounds for volume to be read in and searched
            if sliceLow<1 && periodicZ_==0
                sliceLow=sliceStart;
            end
            sliceHigh=ceil(slice+maxRad);
            if sliceHigh>sliceEnd && periodicZ_==0
                sliceHigh=sliceEnd;
            end
            sizeVolZ=sliceHigh-sliceLow+1;
            
            kkcount=0;
            for kk = sliceLow:sliceHigh
                kkcount = kkcount+1;
                if periodicZ_==1
                    sliceNum=mod(kk-1,lz)+1;
                else
                    sliceNum=kk;
                end
                readVolforSv(:,:,kkcount)=double(0<imread([ tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceNum,'%0.4i')],'tiff'));
            end
            
            sizeX=size(readVolforSv,1);
            sizeY=size(readVolforSv,2);
            sizeZ=size(readVolforSv,3);
            
            if addwallx_==1
                readVolforSv(1,:,:)=ones(1,sizeY,sizeZ);
                readVolforSv(sizeX,:,:)=ones(1,sizeY,sizeZ);
            end
            if addwally_==1
                readVolforSv(:,1,:)=ones(sizeX,1,sizeZ);
                readVolforSv(:,sizeY,:)=ones(sizeX,1,sizeZ);
            end
            if periodicY_==1
                jjcount=0;
                for jj=sizeY-maxRad+1:sizeY
                    jjcount=jjcount+1;
                    jjMod=mod(jj,ly)+1;
                    sideVol1(:,jjcount,:)=readVolforSv(:,jjMod,:);
                end
                jjcount=0;
                for jj=1:maxRad
                    jjcount=jjcount+1;
                    jjMod=mod(jj,ly)+1;
                    sideVol2(:,jjcount,:)=readVolforSv(:,jjMod,:);
                end
                tempVolforSv=cat(2,sideVol1,readVolforSv,sideVol2);
                searchVolforSv=ones(sizeX+2,sizeY+maxRad*2+2,sizeVolZ); %add extra solid layers to search volume, i.e. a 'backstop', if unfamiliar google 'baseball backstop', may not be necessary
                searchVolforSv(2:sizeX+1,2:sizeY+maxRad*2+1,:)=tempVolforSv;
            else
                searchVolforSv=ones(sizeX+2,sizeY+2,sizeVolZ); %add extra solid layers to search volume, may not be necessary
                searchVolforSv(2:sizeX+1,2:sizeY+1,:)=readVolforSv;
            end
            
            %%svsqLoc3D(), function
            [svLocVol(:,:,sliceCount), svTot, sqLocVol(:,:,sliceCount), fluidNTot]=svsqLoc3D(searchVolforSv,mfp,periodicY_,maxRad,a,beta);
            
            rSuga(sliceCount)= 2*C1/((6/pi)^0.5 +C1);%% 0.8909,combination parameter of bounce back and full diffusive, Eq.(21), Landry et al.2016
            rLocVol(:,:,sliceCount)=rSuga(sliceCount)*ones(sizeX+2,sizeY+2,1);
            
            svTotAll=svTotAll+svTot;
            fluidNTotAll=fluidNTotAll+fluidNTot;
        end
    else
        %for standard (global viscosity, conventional navier-stokes) LBM sv
        %and sq are user inputs the following values are superseded by values
        %input on command line, these are basically dummies for the lattice
        %file
        sizeX=lx;
        sizeY=ly;
        sizeZ=lz;
        svLocVol=   ones(lx+2,ly+2,lz);%% without slip
        sqLocVol=sq*ones(lx+2,ly+2,lz);
    end
    
    %determine subvolumes to be sent to each processor
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
        zstart(a)=subz*(a-1)+1;
        zend(a-1)=subz*(a-1);
    end
    
    %writeout fluid node type,local sv, local sq, incoming directions for DiffBB to
    %individual files for each processor
    svSum_Sum=0;
    svCount_Sum=0;
    for rank = 0:(numPx*numPy*numPz)-1
        %disp(rank)
        threadk=floor(rank/(numPx*numPy));
        threadj=floor((rank-(threadk*numPx*numPy))/(numPx));
        threadi=rank-(threadj*numPx)-(threadk*numPx*numPy);
        
        fid = fopen([sim_geom_loc '/' fileOutRoot '/' fileOutRoot '_' num2str(rank)], 'w');
        fprintf(fid, '%d\n', uint16(threadi));
        fprintf(fid, '%d\n', uint16(threadj));
        fprintf(fid, '%d\n', uint16(threadk));
        fprintf(fid, '%d\n', uint16(xend(threadi+1)-xstart(threadi+1)+1));
        fprintf(fid, '%d\n', uint16(yend(threadj+1)-ystart(threadj+1)+1));
        fprintf(fid, '%d\n', uint16(zend(threadk+1)-zstart(threadk+1)+1));
        %Determine which ranks are being communicated with
        [commF,sendRank]=rankMPIcomm26(threadi,threadj,threadk,periodicX_,periodicY_,periodicZ_,numPx,numPy,numPz);
        commFcoll(:,rank+1)=commF';
        sendRankcoll(:,rank+1)=sendRank';
        %Face communication
        for pop=1:6
            fprintf(fid, '%d\n', int16(sendRank(pop)));
        end
        %Edge communication
        for pop=7:18
            fprintf(fid, '%d\n', int16(sendRank(pop)));
        end
        %Corner communication
        for pop=19:26
            fprintf(fid, '%d\n', int16(sendRank(pop)));
        end
        
        svSum(rank+1)=0;
        svCount(rank+1)=0;
        for slice = zstart(threadk+1):zend(threadk+1)
            %outSlice=2*ones(lx+2,ly+2);
            if periodicZ_==1
                if slice==sliceStart
                    readVol(:,:,1)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceEnd,'%0.4i')],'tiff'));
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceStart,'%0.4i')],'tiff'));
                    readVol(:,:,3)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceStart+1,'%0.4i')],'tiff'));
                elseif slice==sliceEnd
                    readVol(:,:,1)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceEnd-1,'%0.4i')],'tiff'));
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceEnd,'%0.4i')],'tiff'));
                    readVol(:,:,3)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceStart,'%0.4i')],'tiff'));
                else
                    readVol(:,:,1)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice-1,'%0.4i')],'tiff'));
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice,'%0.4i')],'tiff'));
                    readVol(:,:,3)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice+1,'%0.4i')],'tiff'));
                end
            else
                if slice==sliceStart
                    readVol(:,:,1)=ones(lx,ly,1);
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceStart,'%0.4i')],'tiff'));
                    readVol(:,:,3)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceStart+1,'%0.4i')],'tiff'));
                elseif slice==sliceEnd
                    readVol(:,:,1)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceEnd-1,'%0.4i')],'tiff'));
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(sliceEnd,'%0.4i')],'tiff'));
                    readVol(:,:,3)=ones(lx,ly,1);
                else
                    readVol(:,:,1)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice-1,'%0.4i')],'tiff'));
                    readVol(:,:,2)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice,'%0.4i')],'tiff'));
                    readVol(:,:,3)=double(0<imread([tiff_dir '/' sample_name '/' sample_name '_' num2str(slice+1,'%0.4i')],'tiff'));
                end
            end
            
            if addwallx_==1
                readVol(1,:,:)=ones(1,ly,3);
                readVol(lx,:,:)=ones(1,ly,3);
            end
            if addwally_==1
                readVol(:,1,:)=ones(lx,1,3);
                readVol(:,ly,:)=ones(lx,1,3);
            end
            
            if periodicY_==1
                searchVol=ones(sizeX+2,sizeY+2,3);
                tempVol=cat(2,readVol(:,ly,:),readVol,readVol(:,1,:));
                searchVol(2:sizeX+1,:,:)=tempVol;
            else
                searchVol=ones(sizeX+2,sizeY+2,3);
                searchVol(2:sizeX+1,2:sizeY+1,:)=readVol;
            end
            
            
            for j=ystart(threadj+1)+1:yend(threadj+1)+1
                for i=xstart(threadi+1)+1:xend(threadi+1)+1
                    sumTemp=0;
                    for pop=1:18
                        sumTemp=sumTemp+searchVol(i+cx(pop),j+cy(pop),2+cz(pop));
                    end
                    if searchVol(i,j,2)==0
                        outSlice(i,j)=0;
                        fprintf(fid, '%d\n', uint8(0));                    %nodetype=0, fluid node
                        fprintf(fid, '%f\n', double(svLocVol(i,j,slice-sliceStart+1))); %sv
                        fprintf(fid, '%f\n', double(sqLocVol(i,j,slice-sliceStart+1))); %sq
                        svCount(rank+1)=svCount(rank+1)+1;
                        svSum(rank+1) = svSum(rank+1) + svLocVol(i,j,slice-sliceStart+1);
                    elseif sumTemp==18
                        outSlice(i,j)=2;
                        fprintf(fid, '%d\n', uint8(2));                    %nodetype=2, solid node
                    else
                        outSlice(i,j)=1;
                        fprintf(fid, '%d\n', uint8(1));                    %nodetype=1, bb node
                        %fprintf(fid, '%f\n', double(rLocVol(i,j,slice-sliceStart+1)));  %r
                        for pop = 1:18                                     %check if BB will receive stream from each pop
                            if searchVol(i-cx(pop),j-cy(pop),2-cz(pop))==0
                                ecomm=1;
                            else
                                ecomm=0;
                            end
                            fprintf(fid, '%d\n' ,uint8(ecomm));
                        end
                    end
                end
            end
            %imagesc(outSlice),drawnow;
        end
        fclose(fid);
        svMean(rank+1)=svSum(rank+1)/svCount(rank+1);
        
        svSum_Sum=svSum_Sum+svSum(rank+1);
        svCount_Sum=svCount_Sum+svCount(rank+1);
        
    end
    svMean=svTotAll/fluidNTotAll;
    svMean_Check=svSum_Sum/svCount_Sum;
    
    
    
    
    %-----To get normalized mfp----------
    MFPLoc=svLocVol(2:lx+1,2:ly+1,:);
    MFPLoc=(1/MFPLoc-0.5)*(pi/6)^0.5;
    MFPLoc=MFPLoc/mfp;
    save ([mfp_loc '/' sample_name '_' num2str(pressure) '_MPa'], 'MFPLoc');
    %--------------------------------------
    
    
    % get slurmjob file
    writeSlurmJob(sample_name,  fileOutRoot , ...
                    pressure, sim_size, sim_geom_loc)
end
end

