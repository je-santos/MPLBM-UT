function [svLoc,svTot,sqLoc,fluidTot] = svsqLoc3D(Vol,mfp,perY,maxRad,alpha,beta)


C1=1.11;
C2=0.61;
        
cx  = [1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0];
cy  = [0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1];
cz  = [0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1];
opp = [2,  1,  4,  3,  6,  5,  10, 9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15];

factX=((cx.^2+cy.^2+cz.^2).^0.5).*cx; %convert radius to cartesian
factY=((cx.^2+cy.^2+cz.^2).^0.5).*cy;
factZ=((cx.^2+cy.^2+cz.^2).^0.5).*cz;

for k = 1:18
    if factX(k)~=0  %% ~= represents !=
        factX(k)=1/factX(k);
    end
    if factY(k)~=0
        factY(k)=1/factY(k);
    end
    if factZ(k)~=0
        factZ(k)=1/factZ(k);
    end
end
sizeX=size(Vol,1);
sizeY=size(Vol,2);
sizeZ=size(Vol,3);
if perY==1
    svLoc=zeros(sizeX,sizeY-maxRad*2);
    sqLoc=zeros(sizeX,sizeY-maxRad*2);
    psiImage=zeros(sizeX,sizeY-maxRad*2);
   psiImage=ones(sizeX,sizeY-maxRad*2);        
else
    svLoc=zeros(sizeX,sizeY);
    sqLoc=zeros(sizeX,sizeY);
    psiImage=zeros(sizeX,sizeY);  
end

k=round(sizeZ/2);

fluidTot=0;
psiTot=0;
if perY==1
    startJ=maxRad+1;
    endJ=sizeY-maxRad;
else
    startJ=1;
    endJ=sizeY;
end
counti=0;
for i = 1:sizeX
    counti=counti+1;
    countj=0;
    for j = startJ:endJ
        countj=countj+1;
        if Vol(i,j,k)==0    %% if Vol is fluid nodes, sv is nonzero, otherwise sv=0
            fluidTot=fluidTot+1;
            %find location of wall 
            %psi = zeros(1,18);
            for discR = 1:18
                check=0;
                r=0;
                while check==0
                    r=r+0.5;
                    jj=j+round(r*factY(discR));
                    ii=i+round(r*factX(discR));
                    kk=k+round(r*factZ(discR));
                    if jj>sizeY || ii>sizeX || kk>sizeZ || kk<1 || jj<1 || ii<1      
                        check=1;
                        r=1000; %dummy, results in a psi=1.000000
                    elseif Vol(ii,jj,kk)==1
                        check=1;
                    end
                end
                d(discR)=r;
            end
            effMFP=0;
            for discR = 1:18
                Href=d(discR)+d(opp(discR));
                PrLarger=exp(-Href/mfp);
                PsiSmaller=(2/pi)*atan(alpha*(2^0.5)*(d(discR)/mfp)^beta); %currently psi Small and large are the same
                PsiLarger=(2/pi)*atan(alpha*(2^0.5)*(d(discR)/mfp)^beta);
                effMFP=effMFP+(1/18)*((1-PrLarger)*PsiSmaller*mfp + PrLarger*PsiLarger*Href);
            end
            svLoc(counti,countj)=1/(0.5+((6/pi)^0.5)*effMFP);       
        end
    end
end

%imagesc(svLoc), drawnow, pause(0.2) 
svTot=sum(sum(svLoc));

%determine local sq using equation 34 of Suga et al. 2011
tauSv=1./svLoc;
tauQtildeSuga=(3+pi*((2*tauSv-1).^2)*C2)./(8*(2*tauSv-1));
sqLocTemp=1./(tauQtildeSuga+0.5);
%sqLocTemp=8.*(2-svLoc)./(8-svLoc);
sqLoc=(sqLocTemp.*(svLoc>0)); %remove infs, doesnt so this need fix 


end

