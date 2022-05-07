
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "mpi.h"



int main(int argc,char* argv[]){

    int size, rank;
    double walltime1, walltime2;
    MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    walltime1=MPI_Wtime();

    double pi = 3.14159265359;
    double cs = pow((1./3.),0.5);
    double sqrt2 = pow(2.,0.5);


///***************BEGIN COMMAND LINE USER INPUTS*****************************************************************************************//
///*************************************************************************************************************************//

    //read simulation name
        char nameFile[150];
    sscanf(argv[1],"%s",nameFile);
    if (rank==0)printf("Simulation Name: %s\n", nameFile);
    //read geometry file name from command line
    char geomFile[150];
    sscanf(argv[2],"%s",geomFile);
	mkdir(argv[2], S_IRWXG);	// creat output folder
    //user input number of processors along each dimension
    const int numPx = atoi(argv[3]);
    const int numPy = atoi(argv[4]);
    const int numPz = atoi(argv[5]);
    if (rank==0)printf("number of processors along x-direction: %d\n", numPx);
    if (rank==0)printf("number of processors along y-direction: %d\n", numPy);
    if (rank==0)printf("number of processors along z-direction: %d\n", numPz);
    if (rank==0)printf("MPICOMM Size: %d\n",size);

    // lattice size
    const int lx = atoi(argv[6]);
    const int ly = atoi(argv[7]);
    const int lz = atoi(argv[8]);

    int flowDir=2;            //direction of flow 0 = x, 1 = y, 2 = z, when no flow default to 2

    //declare periodic boundaries
    int periodic_x = atoi(argv[9]);
    int periodic_y = atoi(argv[10]);
    int periodic_z = atoi(argv[11]);

    double accell[3] = {atof(argv[12]), atof(argv[13]), atof(argv[14])};    //external forces

    //declare pressure boundaries with option of velocity boundaries
    // Zou/He 1997 type using CORRECTED Hecht/Harting 2010 construction
    int pbound_X = atoi(argv[15]);
    int pbound_X_LinP = atoi(argv[16]);
    double pInX = atoi(argv[17]);
    double pOutX = atoi(argv[18]);
    double pInX_uy = 0.0;
    double pInX_uz = 0.0;
    double pOutX_uy = 0.0;
    double pOutX_uz = 0.0;
    double NyX, NzX;

    int pbound_Y = atoi(argv[19]);
    int pbound_Y_LinP = atoi(argv[20]);
    double pInY = atoi(argv[21]);
    double pOutY = atoi(argv[22]);
    double pInY_ux = 0.0;
    double pInY_uz = 0.0;
    double pOutY_ux = 0.0;
    double pOutY_uz = 0.0;
    double NxY, NzY;

    int pbound_Z = atoi(argv[23]);
    int pbound_Z_LinP = atoi(argv[24]);
    double pInZ = atof(argv[25]);
    double pOutZ = atof(argv[26]);
    double pInZ_ux = 0.0;
    double pInZ_uy = 0.0;
    double pOutZ_ux = 0.0;
    double pOutZ_uy = 0.0;
    double NxZ, NyZ;

     /// reservoir boundaries constant density and velocity
    int resbound_X = atoi(argv[27]);
    double resRhoInX = atof(argv[28]);
    double resRhoOutX = atof(argv[29]);

    int resbound_Y = atoi(argv[30]);
    double resRhoInY = atof(argv[31]);
    double resRhoOutY = atof(argv[32]);

    int resbound_Z = atoi(argv[33]);
    double resRhoInZ = atof(argv[34]);
    double resRhoOutZ = atof(argv[35]);

    ///D3Q19
    int sizepop=19;

    ///type of isotropy, enter 1 for original 4th order, 2 for 8th order, 3 for 10th order (untested)
    const int isoNum = atoi(argv[36]);

    /// set globalVisc_ to 1 to turn off local effective viscosity, the kinematic visc is read from command line, code still requires LEV/DiffBB lattice file
    int globalVisc_ = atoi(argv[37]);
    double mu;
    double nu = atof(argv[38]);
    double mfp = atof(argv[39]); ///used to calculate the unbound kinematic viscosity

    ///wall velocity only in Z direction here, but easily changed
    int stdBB_= atoi(argv[40]); ///set to 1 to turn off diffusive bunceback, code still requires LEV/DiffBB lattice file
	double uZwall = atof(argv[41]); /*???   */
    double uWall[3]={0,0,uZwall};	/*   */
	double rBB = atof(argv[42]);	/*rBB=0.8909*/

    ///Shan-Chen SCMP pseudo-potential model
    if(rank==0)printf("Using SC EOS\n");
    int pop, popRS, popM;
    double forceSCx, forceSCy, forceSCz, psiSC1, psiSC2;
    double potGSC,psiC,rhoSC;
    potGSC = atof(argv[43]); //choices based on yu and fan 2010,  GSC='-4.1'
	/*potGSC fuild-fulid?*/
    double rhoTemp, rhoTemp1, rhoTemp2, rhoTempInv, uxTemp, uyTemp, uzTemp, forceAllxTemp, forceAllyTemp, forceAllzTemp;;

    ///defaults to random density, the density entries are superseded if rho file is read
    double rhoFluid = atof(argv[44]);
    double rhoFluidran = atof(argv[45]);
    double rhoWall = atof(argv[46]);
    mu=nu*rhoFluid;
    ///option to read in rho
    int readRhoIn = atoi(argv[47]);
    char rhoFile[150];
    sscanf(argv[48],"%s",rhoFile);
    if(readRhoIn==0){
        if(rank==0) printf("rhoFluid = %lf, rhoFluidran = %lf, rhoWall = %lf, G = %lf\n", rhoFluid,rhoFluidran,rhoWall,potGSC);
    }
    else{
        if(rank==0) printf("G = %lf, read rho file in: %s\n",potGSC,rhoFile);
    }


   ///convergence criteria
    double epsilon = atof(argv[49]);
    double uDarcy_old = 100000;
    double rho_old = 100000;
    char typeStrConv[3]; //only two options now, 'vel' for velocity convergence and 'rho' for density convergence
    sscanf(argv[50],"%s",typeStrConv);
    int typeConv;

    if (strcmp("vel",typeStrConv)==0) typeConv=0;
    else if (strcmp("rho",typeStrConv)==0) typeConv=1;
    else {
            printf("type of convergence '%s' not recognized, default to velocity convergence\n", typeStrConv);
            typeConv=0;
    }

    ///analysis
    double uAll, perm, permFlow, convU, convRho, gradP_perm;
    gradP_perm=1.0; //dummy pressure gradients for now, perm calc different from single phase.
    double uDarcy_mean, uFlow_mean;
    double uDarcy_sum=0;
    double rhoAll, rho_mean;
    double rho_sum=0;
    int thread;

    //Iteration limit and output iters
    int iT, previT;
    int maxT = atoi(argv[51]);
    int tConv = atoi(argv[52]);
    int tOut = atoi(argv[53]);
    double outTemp;
    int outTempInt;
    int countOut;
    int outux_ = atoi(argv[54]);
    int outuy_ = atoi(argv[55]);
    int outuz_ = atoi(argv[56]);
    int outrho_ = atoi(argv[57]);
    int outfIn_ = atoi(argv[58]);
    int outfOut_ = atoi(argv[59]);
    int outgeom_ = atoi(argv[60]);

///********************END USER INPUTS*************************************************************************************************************************//
///************************************************************************************************************************************************************//

///+++++++++++++++++++++Read individual lattice files into each processor++++++++++++++++++++///
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///

        ///read in location of each processor and size of processor lattice subdomain
    int threadi, threadj, threadk;
    int sizeX, sizeY, sizeZ;
    int sizeX2, sizeY2, sizeZ2;
    int commRType;
    int kk;
    char geomFileO[350];

    sprintf(geomFileO, "%s/%s_%d", geomFile, geomFile, rank);
    FILE *geom = fopen(geomFileO,"r");
    if (geom == NULL) {
        perror("error opening geometry file");
    }
    fscanf(geom,"%d",&threadi);
    fscanf(geom,"%d",&threadj);
    fscanf(geom,"%d",&threadk);
    fscanf(geom,"%d",&sizeX);
    fscanf(geom,"%d",&sizeY);
    fscanf(geom,"%d",&sizeZ);
    sizeX2=sizeX+2*isoNum;
    sizeY2=sizeY+2*isoNum;
    sizeZ2=sizeZ+2*isoNum;
    int numNeighbor=18;
    if(isoNum>1){numNeighbor=26;}

    int commRank[26];
    for(kk=0;kk<26;kk++){
        fscanf(geom,"%d",&commRank[kk]); // -1 indicates no communication
    }

    ///input geometry into each processor //3 = ghost node, 2 = interior node, 1 = bb node, 0 = fluid node
    int *bbRegionSub, *indexFunc, *fluidCommSub;
    bbRegionSub=(int *)malloc(sizeof(int)*sizeX2*sizeY2*sizeZ2);
    fluidCommSub=(int *)malloc(sizeof(int)*sizepop*sizeX2*sizeY2*sizeZ2); //indicates if a bb node receives the particle distribution function from any particular velocity e
    indexFunc=(int *)malloc(sizeof(int)*sizeX2*sizeY2*sizeZ2); //indexfunction is 1 for solid and bb, 0 for fluid, unnecessary for now
    double *svLocSub, *sqLocSub;
    svLocSub=(double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    sqLocSub=(double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);


    ///Declare variables for loops
    int k,j,i;
    double indexTemp;
    int tempID, tempID2, tempIDpop, tempIDpop2;

    ///first fill with zeros, 1s, or 3s
    for(k=0;k<sizeZ2;k++){
        for(j=0;j<sizeY2;j++){
            for(i=0;i<sizeX2;i++){
                tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                bbRegionSub[tempID]=3;
                indexFunc[tempID]=1;
                svLocSub[tempID]=0.;
                sqLocSub[tempID]=0.;
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(sizepop*tempID)+pop;
                    fluidCommSub[tempIDpop]=0;
                }
            }
        }
    }
    double svCount, svSum, sqCount, sqSum, svMean, sqMean, svUnbound;
    ///read node type
   for(k=isoNum;k<sizeZ2-isoNum;k++){
        for(j=isoNum;j<sizeY2-isoNum;j++){
            for(i=isoNum;i<sizeX2-isoNum;i++){
                tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                fscanf(geom,"%d",&bbRegionSub[tempID]);
                if(bbRegionSub[tempID]==0){
                    fscanf(geom,"%lf",&svLocSub[tempID]);
                    fscanf(geom,"%lf",&sqLocSub[tempID]);
                    svCount=svCount+1.;
                    svSum = svSum + svLocSub[tempID];
                    sqCount=sqCount+1.;
                    sqSum = sqSum + sqLocSub[tempID];
                }
                if(bbRegionSub[tempID]==1){
                    for(pop=1;pop<=18;pop++){
                        tempIDpop=(sizepop*tempID)+pop;
                        fscanf(geom,"%d",&fluidCommSub[tempIDpop]);
                    }
                }
            }
        }
    }

    fclose(geom);

///+++++++++++++++++++++++++++++++Declare model parameters+++++++++++++++++++++++++++++++///
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///

    ///declare macroscopic variable arrays with ghost nodes
    double *rho, *ux, *uy, *uz;
    rho = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    ux = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    uy = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    uz = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    double uFx,uFy,uFz, Fx,Fy,Fz, MF[19];
    double cxfIn, cyfIn, czfIn, cu, jx, jy, jz;
    double *forceAllx, *forceAlly, *forceAllz;
    forceAllx = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    forceAlly = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    forceAllz = (double *)malloc(sizeof(double)*sizeX2*sizeY2*sizeZ2);
    ///Particle distribution function arrays:
    double *fOut, *fIn;
    fOut = (double *)malloc(sizeof(double)*sizepop*sizeX2*sizeY2*sizeZ2);
    fIn = (double *)malloc(sizeof(double)*sizepop*sizeX2*sizeY2*sizeZ2);

    double An, fSum, sigma, mEq[19], m[19], colM[19], f[19], source[19], colMF[19];
    //D3Q19 velocity model
    double fCalc[19];
    double t[19] = {1./3.,  1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.}; //weights

    int cx[19]  =   {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
    int cy[19]  =   {0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1};
    int cz[19]  =   {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1};

    double cxD[19]  =   {0.,  1., -1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.,  1., -1.,  1., -1.,  0.,  0.,  0.,  0.};
    double cyD[19]  =   {0.,  0.,  0.,  1., -1.,  0.,  0.,  1.,  1., -1., -1.,  0.,  0.,  0.,  0.,  1., -1.,  1., -1.};
    double czD[19]  =   {0.,  0.,  0.,  0.,  0.,  1., -1.,  0.,  0.,  0.,  0.,  1.,  1., -1., -1.,  1.,  1., -1., -1.};


    int opp[19] =   {0,  2,  1,  4,  3,  6,  5,  10, 9,  8,  7, 14, 13, 12, 11, 18, 17, 16, 15};

    double srho, se, sEp, sj, sPi, sm, omEp, omxx, omEpj, sv, sq;
    if(globalVisc_==0){
        ///Guo 2008 parameterization for LEV
        srho = 0.0;
        se = 1.19;
        sEp = 1.4;
        sj = srho;
        sPi = sEp;
        sm = 1.98;
        omEp = 3.0;
        omxx = -0.5;
        omEpj = -11./2.;
        sv = 1.0; //dummy variable local value will be read in
        sq = 1.0; //dummy variable local value will be read in
    }
    else{
        ///Relaxation times following Ammar et al. 2017 used for global viscosity
        sv = 1./(3.*nu+0.5);
        srho = 0.0;
        se = 1.1;
        sEp = 1.1;
        sj = srho;
        sPi = 1.0;
        sq = 1.1;
        sm = sq;
        omEp = 3.;
        omxx = -1./2.;
        omEpj = -11./2.;
    }
    double S[19] = {srho, se, sEp, sj, sq, sj, sq, sj, sq, sv, sPi, sv, sPi, sv, sv, sv, sm, sm, sm};
    //MRT transformation matrices
    double M[19][19] = {
                        {
                             1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0
                        },
                        {
                           -30.0,-11.0,-11.0,-11.0,-11.0,-11.0,-11.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0
                        },
                        {
                            12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0
                        },
                        {
                             0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0
                        },
                        {
                             0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0
                        },
                        {
                             0.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0
                        },
                        {
                             0.0, -4.0, -4.0,  2.0,  2.0,  2.0,  2.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0
                        },
                        {
                             0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0, -2.0, -2.0,  2.0,  2.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0, -1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0
                        },
                        {
                             0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0
                        }
                        };
    double invM[19][19] = {
                        {
                             1./19.,   -5./399.,  1./21.,     0.0,     0.0,      0.0,     0.0,     0.0,     0.0,      0.0,     0.0,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63.,  1./10., -1./10.,      0.0,     0.0,     0.0,     0.0,   1./18., -1./18.,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63., -1./10.,  1./10.,      0.0,     0.0,     0.0,     0.0,   1./18., -1./18.,     0.0,     0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63.,     0.0,     0.0,   1./10., -1./10.,     0.0,     0.0,  -1./36.,  1./36.,  1./12., -1./12.,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63.,     0.0,     0.0,  -1./10.,  1./10.,     0.0,     0.0,  -1./36.,  1./36.,  1./12., -1./12.,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63.,     0.0,     0.0,      0.0,     0.0,  1./10., -1./10.,  -1./36.,  1./36., -1./12.,  1./12.,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19., -11./2394., -1./63.,     0.0,     0.0,      0.0,     0.0, -1./10.,  1./10.,  -1./36.,  1./36., -1./12.,  1./12.,    0.0,    0.0,    0.0,    0.0,    0.0,   0.0
                        },
                        {
                             1./19.,   4./1197., 1./252.,  1./10.,  1./40.,   1./10.,  1./40.,     0.0,     0.0,   1./36.,  1./72.,  1./12.,  1./24.,  1./4.,    0.0,    0.0,  1./8., -1./8.,   0.0
                        },
                        {
                             1./19.,   4./1197., 1./252., -1./10., -1./40.,   1./10.,  1./40.,     0.0,     0.0,   1./36.,  1./72.,  1./12.,  1./24., -1./4.,    0.0,    0.0, -1./8., -1./8.,   0.0
                        },
                        {
                             1./19.,   4./1197., 1./252.,  1./10.,  1./40.,  -1./10., -1./40.,     0.0,     0.0,   1./36.,  1./72.,  1./12.,  1./24., -1./4.,    0.0,    0.0,  1./8.,  1./8.,   0.0
                        },
                        {
                             1./19.,   4./1197., 1./252., -1./10., -1./40.,  -1./10., -1./40.,     0.0,     0.0,   1./36.,  1./72.,  1./12.,  1./24.,  1./4.,    0.0,    0.0, -1./8.,  1./8.,   0.0
                        },
                        {
                             1./19.,   4./1197., 1./252.,  1./10.,  1./40.,      0.0,     0.0,  1./10.,  1./40.,   1./36.,  1./72., -1./12., -1./24.,    0.0,    0.0,  1./4., -1./8.,    0.0,  1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252., -1./10., -1./40.,      0.0,     0.0,  1./10.,  1./40.,   1./36.,  1./72., -1./12., -1./24.,    0.0,    0.0, -1./4.,  1./8.,    0.0,  1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252.,  1./10.,  1./40.,      0.0,     0.0, -1./10., -1./40.,   1./36.,  1./72., -1./12., -1./24.,    0.0,    0.0, -1./4., -1./8.,    0.0, -1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252., -1./10., -1./40.,      0.0,     0.0, -1./10., -1./40.,   1./36.,  1./72., -1./12., -1./24.,    0.0,    0.0,  1./4.,  1./8.,    0.0, -1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252.,     0.0,     0.0,   1./10.,  1./40.,  1./10.,  1./40.,  -1./18., -1./36.,    0.0,      0.0,    0.0,  1./4.,    0.0,    0.0,  1./8., -1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252.,     0.0,     0.0,  -1./10., -1./40.,  1./10.,  1./40.,  -1./18., -1./36.,    0.0,      0.0,    0.0, -1./4.,    0.0,    0.0, -1./8., -1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252.,     0.0,     0.0,   1./10.,  1./40., -1./10., -1./40.,  -1./18., -1./36.,    0.0,      0.0,    0.0, -1./4.,    0.0,    0.0,  1./8.,  1./8.
                        },
                        {
                             1./19.,   4./1197., 1./252.,     0.0,     0.0,  -1./10., -1./40., -1./10., -1./40.,  -1./18., -1./36.,    0.0,      0.0,    0.0,  1./4.,    0.0,    0.0, -1./8.,  1./8.
                        }
                        };

    //wall velocity maxwellian
    //maxwellian equilibrium (same as BGK collision equilibrium fi)
    double cuWall, uWallEq[19];
    for(pop=0;pop<=18;pop++){
        cuWall = cxD[pop]*uWall[0] + cyD[pop]*uWall[1] + czD[pop]*uWall[2];///
        uWallEq[pop]=(1. + 3.*cuWall + (9./2.)*cuWall*cuWall - (3./2.)*(uWall[0]*uWall[0] + uWall[1]*uWall[1] + uWall[2]*uWall[2]));
        if(rank==0)printf("velocity = %i  with uWallEq = %E\n", pop, uWallEq[pop]);
    }

    ///SC order 4 isotropy
    double weightSC[19];
    weightSC[0]=0.;
    for(pop=1;pop<=6;pop++)weightSC[pop]=1./6.; //should total 2 instead of 2/3, some literature describes G*cs*cs*wSC, others G*t, which are equivalent, but former is now more common description
    for(pop=7;pop<=18;pop++)weightSC[pop]=1./12.; //
    rhoSC = 1.0;

    ///SC order 8 isotropy
    //8th order isotropy see sbragaglia et al. 2007 and Ammar et al. 2017
    int cx8[92] = {1,-1,0,0,0,0, 1,-1,1,-1,1,-1,1,-1,0,0,0,0, 1,-1,1,-1,1,-1,1,-1, 2,-2,0,0,0,0, 2,-2,2,-2,2,-2,2,-2,1,1,-1,-1,0,0,0,0,1,1,-1,-1,0,0,0,0, 2,-2,2,-2,2,-2,2,-2,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1, 2,-2,2,-2,2,-2,2,-2,0,0,0,0};
    int cy8[92] = {0,0,-1,1,0,0, 1,1,-1,-1,0,0,0,0,1,-1,1,-1, 1,1,-1,-1,1,1,-1,-1, 0,0,2,-2,0,0, 1,1,-1,-1,0,0,0,0,2,-2,2,-2,2,-2,2,-2,0,0,0,0,1,1,-1,-1, 1,1,-1,-1,1,1,-1,-1,2,-2,2,-2,2,-2,2,-2,1,1,1,1,-1,-1,-1,-1, 2,2,-2,-2,0,0,0,0,2,-2,2,-2};
    int cz8[92] = {0,0,0,0,-1,1, 0,0,0,0,1,1,-1,-1,1,1,-1,-1, 1,1,1,1,-1,-1,-1,-1, 0,0,0,0,2,-2, 0,0,0,0,1,1,-1,-1,0,0,0,0,1,1,-1,-1,2,-2,2,-2,2,-2,2,-2, 1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,2,-2,2,-2,2,-2,2,-2, 0,0,0,0,2,2,-2,-2,2,2,-2,-2};
    double weight8[92];
    for(k=0;k<=5;k++){weight8[k]=4./45.;}
    for(k=6;k<=17;k++){weight8[k]=1./21.;}
    for(k=18;k<=25;k++){weight8[k]=2./105.;}
    for(k=26;k<=31;k++){weight8[k]=5./504.;}
    for(k=32;k<=55;k++){weight8[k]=1./315.;}
    for(k=56;k<=79;k++){weight8[k]=1./630.;}
    for(k=80;k<=91;k++){weight8[k]=1./5040.;}

    ///SC order 10 isotropy
    //10th order isotropy see sbragaglia et al. 2007 and Ammar et al. 2017
    int cx10[170] = {1,-1,0,0,0,0, 1,-1,1,-1,1,-1,1,-1,0,0,0,0, 1,-1,1,-1,1,-1,1,-1, 2,-2,0,0,0,0, 2,-2,2,-2,2,-2,2,-2,1,1,-1,-1,0,0,0,0,1,1,-1,-1,0,0,0,0, 2,-2,2,-2,2,-2,2,-2,1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1, 2,-2,2,-2,2,-2,2,-2,0,0,0,0, 2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,2,-2,1,1,1,1,-1,-1,-1,-1, 3,-3,0,0,0,0, 3,-3,3,-3,3,-3,3,-3,1,1,-1,-1,0,0,0,  0,0,0,0,  0,1,1,-1,-1, 3,-3,3,-3,3,-3,3,-3,1,1,-1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1};
    int cy10[170] = {0,0,-1,1,0,0, 1,1,-1,-1,0,0,0,0,1,-1,1,-1, 1,1,-1,-1,1,1,-1,-1, 0,0,2,-2,0,0, 1,1,-1,-1,0,0,0,0,2,-2,2,-2,2,-2,2,-2,0,0,0,0,1,1,-1,-1, 1,1,-1,-1,1,1,-1,-1,2,-2,2,-2,2,-2,2,-2,1,1,-1,-1,1,1,-1,-1, 2,2,-2,-2,0,0,0,0,2,-2,2,-2, 2,2,-2,-2,2,2,-2,-2,1,1,1,1,-1,-1,-1,-1,2,-2,2,-2,2,-2,2,-2, 0,0,3,-3,0,0, 1,1,-1,-1,0,0,0,  0,3,-3,3,-3,3,-3,3,-3,1,1,-1,-1,0,0,0,  0, 1,1,-1,-1,1,1,-1,-1,3,-3,3,-3,3,-3,3,-3,1,1,-1,-1,1,1,-1,-1};
    int cz10[170] = {0,0,0,0,-1,1, 0,0,0,0,1,1,-1,-1,1,1,-1,-1, 1,1,1,1,-1,-1,-1,-1, 0,0,0,0,2,-2, 0,0,0,0,1,1,-1,-1,0,0,0,0,1,1,-1,-1,2,-2,2,-2,2,-2,2,-2, 1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,2,-2,2,-2,2,-2,2,-2, 0,0,0,0,2,2,-2,-2,2,2,-2,-2, 1,1,1,1,-1,-1,-1,-1,2,2,-2,-2,2,2,-2,-2,2,2,-2,-2,2,2,-2,-2, 0,0,0,0,3,-3, 0,0,0,  0,1,1,-1,-1,0,0,0,  0,1,1,-1,-1,3,-3,3,-3,3,-3,3,-3, 1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,3,-3,3,-3,3,-3,3,-3};
    double weight10[170];
    for(k=0;k<=5;k++){weight10[k]=353./5355.;}
    for(k=6;k<=17;k++){weight10[k]=38./1071.;}
    for(k=18;k<=25;k++){weight10[k]=271./14280.;}
    for(k=26;k<=31;k++){weight10[k]=139./14280.;}
    for(k=32;k<=55;k++){weight10[k]=53./10710.;}
    for(k=56;k<=79;k++){weight10[k]=5./2142.;}
    for(k=80;k<=91;k++){weight10[k]=41./85680.;}
    for(k=92;k<=115;k++){weight10[k]=1./4284.;}
    for(k=116;k<=121;k++){weight10[k]=1./5355.;}
    for(k=122;k<=145;k++){weight10[k]=1./10710.;}
    for(k=146;k<=169;k++){weight10[k]=1./42840.;}


///+++++++++++++++++++++++++++Declare variables for D3Q19 stream/scalar MPI passing+++++++++++++++++++++++++///
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///
    double *fIn_send, *fIn_recv, *scalar_send, *scalar_recv;
    int *ind_send, *ind_recv;
    fIn_send = (double *)malloc(sizeof(double)*sizeX*sizeY*2*5);
    fIn_recv = (double *)malloc(sizeof(double)*sizeX*sizeY*2*5);
    scalar_send = (double *)malloc(sizeof(double)*sizeX*sizeY*(1+isoNum));
    scalar_recv = (double *)malloc(sizeof(double)*sizeX*sizeY*(1+isoNum));
    ind_send = (int *)malloc(sizeof(int)*sizeX*sizeY*(1+isoNum));
    ind_recv = (int *)malloc(sizeof(int)*sizeX*sizeY*(1+isoNum));

    int oppComm[26] = {1, 0, 3, 2, 5, 4, 9, 8, 7, 6, 13, 12, 11, 10, 17, 16, 15, 14, 25, 24, 23, 22, 21, 20, 19, 18}; //commtype = pop-1, only use comms 0-17 for streaming, corner pass for rho pass
    int countS, tagS, tagR, tag, eDir;
    int passDir[18], sizeSend[18];
    int numComm=0;
    int numPop[18] = {5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,1,1};
    int numCommR=0;
    int passDirR[numNeighbor], sizeRSend[numNeighbor];

    int passType[26]; //determine type of complimentary pass, either eDir = + then - OR - then +, depending on location in subdomain grid
    if(threadi%2==0){passType[0]=0; passType[1]=1;} //face passes
    if(threadi%2==1){passType[0]=1; passType[1]=0;}
    if(threadj%2==0){passType[2]=2; passType[3]=3;}
    if(threadj%2==1){passType[2]=3; passType[3]=2;}
    if(threadk%2==0){passType[4]=4; passType[5]=5;}
    if(threadk%2==1){passType[4]=5; passType[5]=4;}
    if(threadi%2==0){passType[6]=6; passType[7]=9; passType[8]=7; passType[9]=8;} //edge passes
    if(threadi%2==1){passType[6]=9; passType[7]=6; passType[8]=8; passType[9]=7;}
    if(threadk%2==0){passType[10]=10; passType[11]=13; passType[12]=11; passType[13]=12;}
    if(threadk%2==1){passType[10]=13; passType[11]=10; passType[12]=12; passType[13]=11;}
    if(threadj%2==0){passType[14]=14; passType[15]=17; passType[16]=15; passType[17]=16;}
    if(threadj%2==1){passType[14]=17; passType[15]=14; passType[16]=16; passType[17]=15;}
    if(threadk%2==0){passType[18]=18; passType[19]=25; passType[20]=19; passType[21]=24; //corner passes, only used for isoNum>1
                     passType[22]=20; passType[23]=23; passType[24]=21; passType[25]=22;}
    if(threadk%2==1){passType[18]=25; passType[19]=18; passType[20]=24; passType[21]=19;
                     passType[22]=23; passType[23]=20; passType[24]=22; passType[25]=21;}

    ///For passing streamed particle distributions from ghost layer to edge/face of neighboring subdomain D3Q19
    int passPops[18][5] =
                            {
                                    {
                                        1, 7, 9, 11, 13  //+x
                                    },
                                    {
                                        2, 8, 10, 12, 14  //-x
                                    },
                                    {
                                        3, 7, 8, 15, 17  //+y
                                    },
                                    {
                                        4, 9, 10, 16, 18  //-y
                                    },
                                    {
                                        5, 11, 12, 15, 16  //+z
                                    },
                                    {
                                        6, 13, 14, 17, 18  //-z
                                    },
                                    {
                                        7, 0, 0, 0, 0  // e7 (0s are dummy)
                                    },
                                    {
                                        8, 0, 0, 0, 0  // e8
                                    },
                                    {
                                        9, 0, 0, 0, 0  // e9
                                    },
                                    {
                                        10, 0, 0, 0, 0  // e10
                                    },
                                    {
                                        11, 0, 0, 0, 0  // e11
                                    },
                                    {
                                        12, 0, 0, 0, 0  // e12
                                    },
                                    {
                                        13, 0, 0, 0, 0  // e13
                                    },
                                    {
                                        14, 0, 0, 0, 0  // e14
                                    },
                                    {
                                        15, 0, 0, 0, 0  // e15
                                    },
                                    {
                                        16, 0, 0, 0, 0  // e16
                                    },
                                    {
                                        17, 0, 0, 0, 0  // e17
                                    },
                                    {
                                        18, 0, 0, 0, 0  // e18
                                    }
                                };
    int xSendBounds[18][2] =
                            {
                                {
                                   sizeX2-isoNum,sizeX2-isoNum   //+x
                                },
                                {
                                   isoNum-1,isoNum-1   //-x
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //+y
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //-y
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //+z
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //-z
                                },
                                {
                                   sizeX2-isoNum,sizeX2-isoNum   //+x,+y,7
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,+y,8
                                },
                                {
                                   sizeX2-isoNum,sizeX2-isoNum   //+x,-y,9
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,-y,10
                                },
                                {
                                   sizeX2-isoNum,sizeX2-isoNum   //+x,+z,11
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,+z,12
                                },
                                {
                                   sizeX2-isoNum,sizeX2-isoNum   //+x,-z,13
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,-z,14
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //+x,+z,15
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //-x,+z,16
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //+x,-z,17
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //-x,-z,18
                                },
                            };
    int xRecvBounds[18][2] =
                            {
                                {
                                   isoNum,isoNum   //+x
                                },
                                {
                                   sizeX2-1-isoNum,sizeX2-1-isoNum   //-x
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //+y
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //-y
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //+z
                                },
                                {
                                   isoNum,sizeX2-1-isoNum   //-z
                                },
                                {
                                   isoNum,isoNum   //+x,+y,7
                                },
                                {
                                   sizeX2-1-isoNum,sizeX2-1-isoNum   //-x,+y,8
                                },
                                {
                                   isoNum,isoNum   //+x,-y,9
                                },
                                {
                                   sizeX2-1-isoNum,sizeX2-1-isoNum   //-x,-y,10
                                },
                                {
                                   isoNum,isoNum   //+x,+z,11
                                },
                                {
                                   sizeX2-1-isoNum,sizeX2-1-isoNum   //-x,+z,12
                                },
                                {
                                   isoNum,isoNum   //+x,-z,13
                                },
                                {
                                   sizeX2-1-isoNum,sizeX2-1-isoNum   //-x,-z,14
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //+x,+z,15
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //-x,+z,16
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //+x,-z,17
                                },
                                {
                                   isoNum,sizeX2-1-isoNum      //-x,-z,18
                                },
                            };
    int ySendBounds[18][2] =
                            {
                                {
                                   isoNum,sizeY2-1-isoNum   //+x
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x
                                },
                                {
                                   sizeY2-isoNum,sizeY2-isoNum   //+y
                                },
                                {
                                   isoNum-1,isoNum-1   //-y
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+z
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-z
                                },
                                {
                                   sizeY2-isoNum,sizeY2-isoNum   //+x,+y,7
                                },
                                {
                                   sizeY2-isoNum,sizeY2-isoNum   //-x,+y,8
                                },
                                {
                                   isoNum-1,isoNum-1  //+x,-y,9
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,-y,10
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+x,+z,11
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x,+z,12
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+x,-z,13
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x,-z,14
                                },
                                {
                                   sizeY2-isoNum,sizeY2-isoNum      //+x,+z,15
                                },
                                {
                                   isoNum-1,isoNum-1      //-x,+z,16
                                },
                                {
                                   sizeY2-isoNum,sizeY2-isoNum    //+x,-z,17
                                },
                                {
                                   isoNum-1,isoNum-1     //-x,-z,18
                                },
                            };
    int yRecvBounds[18][2] =
                            {
                                {
                                   isoNum,sizeY2-1-isoNum   //+x
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x
                                },
                                {
                                   isoNum,isoNum   //+y
                                },
                                {
                                   sizeY2-1-isoNum,sizeY2-1-isoNum   //-y
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+z
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-z
                                },
                                {
                                   isoNum,isoNum   //+x,+y,7
                                },
                                {
                                   isoNum,isoNum   //-x,+y,8
                                },
                                {
                                   sizeY2-1-isoNum,sizeY2-1-isoNum  //+x,-y,9
                                },
                                {
                                   sizeY2-1-isoNum,sizeY2-1-isoNum   //-x,-y,10
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+x,+z,11
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x,+z,12
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //+x,-z,13
                                },
                                {
                                   isoNum,sizeY2-1-isoNum   //-x,-z,14
                                },
                                {
                                   isoNum,isoNum      //+x,+z,15
                                },
                                {
                                   sizeY2-1-isoNum,sizeY2-1-isoNum      //-x,+z,16
                                },
                                {
                                   isoNum,isoNum    //+x,-z,17
                                },
                                {
                                   sizeY2-1-isoNum,sizeY2-1-isoNum     //-x,-z,18
                                },
                            };
    int zSendBounds[18][2] =
                            {
                                {
                                   isoNum,sizeZ2-1-isoNum   //+x
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //+y
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-y
                                },
                                {
                                   sizeZ2-isoNum,sizeZ2-isoNum   //+z
                                },
                                {
                                   isoNum-1,isoNum-1   //-z
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //+x,+y,7
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x,+y,8
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum  //+x,-y,9
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x,-y,10
                                },
                                {
                                   sizeZ2-isoNum,sizeZ2-isoNum   //+x,+z,11
                                },
                                {
                                   sizeZ2-isoNum,sizeZ2-isoNum   //-x,+z,12
                                },
                                {
                                   isoNum-1,isoNum-1   //+x,-z,13
                                },
                                {
                                   isoNum-1,isoNum-1   //-x,-z,14
                                },
                                {
                                   sizeZ2-isoNum,sizeZ2-isoNum      //+x,+z,15
                                },
                                {
                                   sizeZ2-isoNum,sizeZ2-isoNum      //-x,+z,16
                                },
                                {
                                   isoNum-1,isoNum-1    //+x,-z,17
                                },
                                {
                                   isoNum-1,isoNum-1     //-x,-z,18
                                },
                            };
    int zRecvBounds[18][2] =
                            {
                                {
                                   isoNum,sizeZ2-1-isoNum   //+x
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //+y
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-y
                                },
                                {
                                   isoNum,isoNum   //+z
                                },
                                {
                                   sizeZ2-1-isoNum,sizeZ2-1-isoNum   //-z
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //+x,+y,7
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x,+y,8
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum  //+x,-y,9
                                },
                                {
                                   isoNum,sizeZ2-1-isoNum   //-x,-y,10
                                },
                                {
                                   isoNum,isoNum   //+x,+z,11
                                },
                                {
                                   isoNum,isoNum   //-x,+z,12
                                },
                                {
                                   sizeZ2-1-isoNum,sizeZ2-1-isoNum   //+x,-z,13
                                },
                                {
                                   sizeZ2-1-isoNum,sizeZ2-1-isoNum   //-x,-z,14
                                },
                                {
                                   isoNum,isoNum      //+x,+z,15
                                },
                                {
                                   isoNum,isoNum      //-x,+z,16
                                },
                                {
                                   sizeZ2-1-isoNum,sizeZ2-1-isoNum    //+x,-z,17
                                },
                                {
                                   sizeZ2-1-isoNum,sizeZ2-1-isoNum     //-x,-z,18
                                },
                            };

    for(kk=0;kk<=17;kk++){
        sizeSend[kk]=numPop[kk]*(zSendBounds[kk][1]-zSendBounds[kk][0]+1)*(ySendBounds[kk][1]-ySendBounds[kk][0]+1)*(xSendBounds[kk][1]-xSendBounds[kk][0]+1);
    }

    ///Create passing sequence D3Q19
    for(kk=0;kk<=17;kk++){
        eDir=passType[kk];
        if (commRank[eDir] > -1){
            numComm=numComm+1;
            passDir[numComm-1]=eDir;
        }
    }

    int boundS[14], boundR[14];
    if(isoNum==1){
        boundS[0] = sizeX; //send inds
        boundS[1] = sizeX;
        boundS[2] = 1;
        boundS[3] = 1;
        boundS[4] = 1;
        boundS[5] = sizeX;

        boundS[6] = sizeY;
        boundS[7] = sizeY;
        boundS[8] = 1;
        boundS[9] = sizeY;

        boundS[10] = sizeZ;
        boundS[11] = sizeZ;
        boundS[12] = 1;
        boundS[13] = sizeZ;

        boundR[0] = 0; //recv inds
        boundR[1] = 0;
        boundR[2] = sizeX+1;
        boundR[3] = sizeX+1;
        boundR[4] = 1;
        boundR[5] = sizeX;

        boundR[6] = sizeY+1;
        boundR[7] = sizeY+1;
        boundR[8] = 1;
        boundR[9] = sizeY;

        boundR[10] = sizeZ+1;
        boundR[11] = sizeZ+1;
        boundR[12] = 1;
        boundR[13] = sizeZ;
    }
    else if(isoNum==2) {
        boundS[0] = sizeX; //send inds
        boundS[1] = sizeX+1;
        boundS[2] = 2;
        boundS[3] = 3;
        boundS[4] = 2;
        boundS[5] = sizeX+1;

        boundS[6] = sizeY;
        boundS[7] = sizeY+1;
        boundS[8] = 2;
        boundS[9] = sizeY+1;

        boundS[10] = sizeZ;
        boundS[11] = sizeZ+1;
        boundS[12] = 2;
        boundS[13] = sizeZ+1;

        boundR[0] = 0; //recv inds
        boundR[1] = 1;
        boundR[2] = sizeX+2;
        boundR[3] = sizeX+3;
        boundR[4] = 2;
        boundR[5] = sizeX+1;

        boundR[6] = sizeY+2;
        boundR[7] = sizeY+3;
        boundR[8] = 2;
        boundR[9] = sizeY+1;

        boundR[10] = sizeZ+2;
        boundR[11] = sizeZ+3;
        boundR[12] = 2;
        boundR[13] = sizeZ+1;
    }
    else {
        boundS[0] = sizeX; //send inds
        boundS[1] = sizeX+2;
        boundS[2] = 3;
        boundS[3] = 5;
        boundS[4] = 3;
        boundS[5] = sizeX+2;

        boundS[6] = sizeY;
        boundS[7] = sizeY+2;
        boundS[8] = 3;
        boundS[9] = sizeY+2;

        boundS[10] = sizeZ;
        boundS[11] = sizeZ+2;
        boundS[12] = 3;
        boundS[13] = sizeZ+2;

        boundR[0] = 0; //recv inds
        boundR[1] = 2;
        boundR[2] = sizeX+3;
        boundR[3] = sizeX+5;
        boundR[4] = 2;
        boundR[5] = sizeX+2;

        boundR[6] = sizeY+3;
        boundR[7] = sizeY+5;
        boundR[8] = 2;
        boundR[9] = sizeY+2;

        boundR[10] = sizeZ+3;
        boundR[11] = sizeZ+5;
        boundR[12] = 2;
        boundR[13] = sizeZ+2;
    }

    ///For passing density/scalar information from corner/edge/face of subdomain to ghost layer of neighboring subdomain
    int xRSendBounds[26][2] =
                        {
                            {
                               boundS[0],boundS[1]   //+x
                            },
                            {
                               boundS[2],boundS[3]   //-x
                            },
                            {
                               boundS[4],boundS[5]   //+y
                            },
                            {
                               boundS[4],boundS[5]   //-y
                            },
                            {
                               boundS[4],boundS[5]   //+z
                            },
                            {
                               boundS[4],boundS[5]   //-z
                            },
                            {
                               boundS[0],boundS[1]   //+x,+y,7
                            },
                            {
                               boundS[2],boundS[3]   //-x,+y,8
                            },
                            {
                               boundS[0],boundS[1]   //+x,-y,9
                            },
                            {
                               boundS[2],boundS[3]   //-x,-y,10
                            },
                            {
                               boundS[0],boundS[1]   //+x,+z,11
                            },
                            {
                               boundS[2],boundS[3]   //-x,+z,12
                            },
                            {
                               boundS[0],boundS[1]   //+x,-z,13
                            },
                            {
                               boundS[2],boundS[3]   //-x,-z,14
                            },
                            {
                               boundS[4],boundS[5]      //+y,+z,15
                            },
                            {
                               boundS[4],boundS[5]      //-y,+z,16
                            },
                            {
                               boundS[4],boundS[5]      //+y,-z,17
                            },
                            {
                               boundS[4],boundS[5]      //-y,-z,18
                            },
                            {
                               boundS[0],boundS[1]      //+x,+y,+z,19
                            },
                            {
                               boundS[2],boundS[3]      //-x,+y,+z,20
                            },
                            {
                               boundS[0],boundS[1]      //x,-y,+z,21
                            },
                            {
                               boundS[2],boundS[3]      //-x,-y,+z,22
                            },
                            {
                               boundS[0],boundS[1]      //+x,+y,-z,23
                            },
                            {
                               boundS[2],boundS[3]      //-x,+y,-z,24
                            },
                            {
                               boundS[0],boundS[1]      //x,-y,-z,25
                            },
                            {
                               boundS[2],boundS[3]      //-x,-y,-z,26
                            }
                        };
    int yRSendBounds[26][2] =
                        {
                            {
                               boundS[8],boundS[9]   //+x
                            },
                            {
                               boundS[8],boundS[9]   //-x
                            },
                            {
                               boundS[6],boundS[7]   //+y
                            },
                            {
                               boundS[2],boundS[3]   //-y
                            },
                            {
                               boundS[8],boundS[9]   //+z
                            },
                            {
                               boundS[8],boundS[9]   //-z
                            },
                            {
                               boundS[6],boundS[7]   //+x,+y,7
                            },
                            {
                               boundS[6],boundS[7]   //-x,+y,8
                            },
                            {
                               boundS[2],boundS[3]  //+x,-y,9
                            },
                            {
                               boundS[2],boundS[3]   //-x,-y,10
                            },
                            {
                               boundS[8],boundS[9]   //+x,+z,11
                            },
                            {
                               boundS[8],boundS[9]   //-x,+z,12
                            },
                            {
                               boundS[8],boundS[9]   //+x,-z,13
                            },
                            {
                               boundS[8],boundS[9]   //-x,-z,14
                            },
                            {
                               boundS[6],boundS[7]      //+y,+z,15
                            },
                            {
                               boundS[2],boundS[3]      //-y,+z,16
                            },
                            {
                               boundS[6],boundS[7]    //+y,-z,17
                            },
                            {
                               boundS[2],boundS[3]     //-y,-z,18
                            },
                            {
                               boundS[6],boundS[7]      //+x,+y,+z,19
                            },
                            {
                               boundS[6],boundS[7]      //-x,+y,+z,20
                            },
                            {
                               boundS[2],boundS[3]      //+x,-y,+z,21
                            },
                            {
                               boundS[2],boundS[3]      //-x,-y,+z,22
                            },
                            {
                               boundS[6],boundS[7]      //+x,+y,-z,23
                            },
                            {
                               boundS[6],boundS[7]      //-x,+y,-z,24
                            },
                            {
                               boundS[2],boundS[3]      //x,-y,-z,25
                            },
                            {
                               boundS[2],boundS[3]      //-x,-y,-z,26
                            }
                        };
    int zRSendBounds[26][2] =
                        {
                            {
                               boundS[12],boundS[13]   //+x
                            },
                            {
                               boundS[12],boundS[13]   //-x
                            },
                            {
                               boundS[12],boundS[13]   //+y
                            },
                            {
                               boundS[12],boundS[13]   //-y
                            },
                            {
                               boundS[10],boundS[11]   //+z
                            },
                            {
                               boundS[2],boundS[3]   //-z
                            },
                            {
                               boundS[12],boundS[13]   //+x,+y,7
                            },
                            {
                               boundS[12],boundS[13]   //-x,+y,8
                            },
                            {
                               boundS[12],boundS[13]  //+x,-y,9
                            },
                            {
                               boundS[12],boundS[13]   //-x,-y,10
                            },
                            {
                               boundS[10],boundS[11]   //+x,+z,11
                            },
                            {
                               boundS[10],boundS[11]   //-x,+z,12
                            },
                            {
                               boundS[2],boundS[3]   //+x,-z,13
                            },
                            {
                               boundS[2],boundS[3]   //-x,-z,14
                            },
                            {
                               boundS[10],boundS[11]      //+y,+z,15
                            },
                            {
                               boundS[10],boundS[11]      //-y,+z,16
                            },
                            {
                               boundS[2],boundS[3]    //+y,-z,17
                            },
                            {
                               boundS[2],boundS[3]     //-y,-z,18
                            },
                            {
                               boundS[10],boundS[11]      //+x,+y,+z,19
                            },
                            {
                               boundS[10],boundS[11]      //-x,+y,+z,20
                            },
                            {
                               boundS[10],boundS[11]      //+x,-y,+z,21
                            },
                            {
                               boundS[10],boundS[11]      //-x,-y,+z,22
                            },
                            {
                               boundS[2],boundS[3]      //+x,+y,-z,23
                            },
                            {
                               boundS[2],boundS[3]      //-x,+y,-z,24
                            },
                            {
                               boundS[2],boundS[3]      //x,-y,-z,25
                            },
                            {
                               boundS[2],boundS[3]      //-x,-y,-z,26
                            }
                        };
    int xRRecvBounds[26][2] =
                        {
                            {
                               boundR[0],boundR[1]   //+x
                            },
                            {
                               boundR[2],boundR[3]   //-x
                            },
                            {
                               boundR[4],boundR[5]   //+y
                            },
                            {
                               boundR[4],boundR[5]   //-y
                            },
                            {
                               boundR[4],boundR[5]   //+z
                            },
                            {
                               boundR[4],boundR[5]   //-z
                            },
                            {
                               boundR[0],boundR[1]   //+x,+y,7
                            },
                            {
                               boundR[2],boundR[3]   //-x,+y,8
                            },
                            {
                               boundR[0],boundR[1]   //+x,-y,9
                            },
                            {
                               boundR[2],boundR[3]   //-x,-y,10
                            },
                            {
                               boundR[0],boundR[1]   //+x,+z,11
                            },
                            {
                               boundR[2],boundR[3]   //-x,+z,12
                            },
                            {
                               boundR[0],boundR[1]   //+x,-z,13
                            },
                            {
                               boundR[2],boundR[3]   //-x,-z,14
                            },
                            {
                               boundR[4],boundR[5]      //+y,+z,15
                            },
                            {
                               boundR[4],boundR[5]      //-y,+z,16
                            },
                            {
                               boundR[4],boundR[5]      //+y,-z,17
                            },
                            {
                               boundR[4],boundR[5]      //-y,-z,18
                            },
                            {
                               boundR[0],boundR[1]      //+x,+y,+z,19
                            },
                            {
                               boundR[2],boundR[3]      //-x,+y,+z,20
                            },
                            {
                               boundR[0],boundR[1]      //+x,-y,+z,21
                            },
                            {
                               boundR[2],boundR[3]      //-x,-y,+z,22
                            },
                            {
                               boundR[0],boundR[1]      //+x,+y,-z,23
                            },
                            {
                               boundR[2],boundR[3]      //-x,+y,-z,24
                            },
                            {
                               boundR[0],boundR[1]      //x,-y,-z,25
                            },
                            {
                               boundR[2],boundR[3]      //-x,-y,-z,26
                            }
                        };
    int yRRecvBounds[26][2] =
                        {
                            {
                               boundR[8],boundR[9]   //+x
                            },
                            {
                               boundR[8],boundR[9]   //-x
                            },
                            {
                               boundR[0],boundR[1]   //+y
                            },
                            {
                               boundR[6],boundR[7]   //-y
                            },
                            {
                               boundR[8],boundR[9]   //+z
                            },
                            {
                               boundR[8],boundR[9]   //-z
                            },
                            {
                               boundR[0],boundR[1]   //+x,+y,7
                            },
                            {
                               boundR[0],boundR[1]   //-x,+y,8
                            },
                            {
                               boundR[6],boundR[7]  //+x,-y,9
                            },
                            {
                               boundR[6],boundR[7]   //-x,-y,10
                            },
                            {
                               boundR[8],boundR[9]   //+x,+z,11
                            },
                            {
                               boundR[8],boundR[9]   //-x,+z,12
                            },
                            {
                               boundR[8],boundR[9]   //+x,-z,13
                            },
                            {
                               boundR[8],boundR[9]   //-x,-z,14
                            },
                            {
                               boundR[0],boundR[1]      //+y,+z,15
                            },
                            {
                               boundR[6],boundR[7]      //-y,+z,16
                            },
                            {
                               boundR[0],boundR[1]    //+y,-z,17
                            },
                            {
                               boundR[6],boundR[7]     //-y,-z,18
                            },
                            {
                               boundR[0],boundR[1]      //+x,+y,+z,19
                            },
                            {
                               boundR[0],boundR[1]      //-x,+y,+z,20
                            },
                            {
                               boundR[6],boundR[7]      //+x,-y,+z,21
                            },
                            {
                               boundR[6],boundR[7]      //-x,-y,+z,22
                            },
                            {
                               boundR[0],boundR[1]      //+x,+y,-z,23
                            },
                            {
                               boundR[0],boundR[1]      //-x,+y,-z,24
                            },
                            {
                               boundR[6],boundR[7]      //x,-y,-z,25
                            },
                            {
                               boundR[6],boundR[7]      //-x,-y,-z,26
                            }
                        };
    int zRRecvBounds[26][2] =
                        {
                            {
                               boundR[12],boundR[13]   //+x
                            },
                            {
                               boundR[12],boundR[13]   //-x
                            },
                            {
                               boundR[12],boundR[13]   //+y
                            },
                            {
                               boundR[12],boundR[13]   //-y
                            },
                            {
                               boundR[0],boundR[1]   //+z
                            },
                            {
                               boundR[10],boundR[11]   //-z
                            },
                            {
                               boundR[12],boundR[13]   //+x,+y,7
                            },
                            {
                               boundR[12],boundR[13]   //-x,+y,8
                            },
                            {
                               boundR[12],boundR[13]  //+x,-y,9
                            },
                            {
                               boundR[12],boundR[13]   //-x,-y,10
                            },
                            {
                               boundR[0],boundR[1]   //+x,+z,11
                            },
                            {
                               boundR[0],boundR[1]   //-x,+z,12
                            },
                            {
                               boundR[10],boundR[11]   //+x,-z,13
                            },
                            {
                               boundR[10],boundR[11]   //-x,-z,14
                            },
                            {
                               boundR[0],boundR[1]      //+y,+z,15
                            },
                            {
                               boundR[0],boundR[1]      //-y,+z,16
                            },
                            {
                               boundR[10],boundR[11]    //+y,-z,17
                            },
                            {
                               boundR[10],boundR[11]     //-y,-z,18
                            },
                             {
                               boundR[0],boundR[1]      //+x,+y,+z,19
                            },
                            {
                               boundR[0],boundR[1]      //-x,+y,+z,20
                            },
                            {
                               boundR[0],boundR[1]      //+x,-y,+z,21
                            },
                            {
                               boundR[0],boundR[1]      //-x,-y,+z,22
                            },
                            {
                               boundR[10],boundR[11]      //+x,+y,-z,23
                            },
                            {
                               boundR[10],boundR[11]      //-x,+y,-z,24
                            },
                            {
                               boundR[10],boundR[11]      //x,-y,-z,25
                            },
                            {
                               boundR[10],boundR[11]      //-x,-y,-z,26
                            }
                        };
    ///size sends for rho/scalar pass
    for(kk=0;kk<numNeighbor;kk++){
        sizeRSend[kk]=(zRSendBounds[kk][1]-zRSendBounds[kk][0]+1)*(yRSendBounds[kk][1]-yRSendBounds[kk][0]+1)*(xRSendBounds[kk][1]-xRSendBounds[kk][0]+1);
    };

    ///Create passing sequence
    for(kk=0;kk<numNeighbor;kk++){
        eDir=passType[kk];
        if (commRank[eDir] > -1){
            numCommR=numCommR+1;
            passDirR[numCommR-1]=eDir;
        }
    }
    if(rank==0){printf("rank 0, isoNum = %d, numCommR = %d\n", isoNum, numCommR);}

///+++++++++++++++++++++Create node type location lists+++++++++++++++++++++++++++++++++///
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///


        ///declare lists to be used in place of conditional rastering
    int *listBBnodes, *listFnodes, *listFintnodes, *listBBPnodes;
    listBBnodes = (int *)malloc(sizeof(int)*sizeX*sizeY*sizeZ);
    listFnodes = (int *)malloc(sizeof(int)*sizeX*sizeY*sizeZ);
    listFintnodes = (int *)malloc(sizeof(int)*sizeX*sizeY*sizeZ);
    listBBPnodes = (int *)malloc(sizeof(int)*sizeX*sizeY*sizeZ);
    int *listPxInnodes, *listPxOutnodes, *listPyInnodes, *listPyOutnodes, *listPzInnodes, *listPzOutnodes;
    listPxInnodes = (int *)malloc(sizeof(int)*sizeY*sizeZ);
    listPxOutnodes = (int *)malloc(sizeof(int)*sizeY*sizeZ);
    listPyInnodes = (int *)malloc(sizeof(int)*sizeX*sizeZ);
    listPyOutnodes = (int *)malloc(sizeof(int)*sizeX*sizeZ);
    listPzInnodes = (int *)malloc(sizeof(int)*sizeY*sizeX);
    listPzOutnodes = (int *)malloc(sizeof(int)*sizeY*sizeX);

    int countPxIn=-1;
    int countPxOut=-1;
    int countPyIn=-1;
    int countPyOut=-1;
    int countPzIn=-1;
    int countPzOut=-1;
    int countBBP=-1;
    int countF=-1;
    int countFint=-1;
    int countBB=-1;
    int iD;

    ///create lists of node locations to be used by main loop, this avoids conditionals in the main loop, note pressure boundary nodes are also fluid nodes
    for(k=isoNum;k<sizeZ2-isoNum;k++){
        for(j=isoNum;j<sizeY2-isoNum;j++){
            for(i=isoNum;i<sizeX2-isoNum;i++){
                tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                if (bbRegionSub[tempID]==0){
                    countF = countF+1;
                    listFnodes[countF]=tempID;
                    if((pbound_X==1 || resbound_X==1) && i==isoNum && threadi==0){
                        countPxIn=countPxIn+1;
                        listPxInnodes[countPxIn]=tempID;
                    }
                    else if((pbound_X==1 || resbound_X==1) && i==sizeX2-1-isoNum && threadi==numPx-1){
                        countPxOut=countPxOut+1;
                        listPxOutnodes[countPxOut]=tempID;
                    }
                    else if((pbound_Y==1 || resbound_Y==1) && j==isoNum && threadj==0){
                        countPyIn=countPyIn+1;
                        listPyInnodes[countPyIn]=tempID;
                    }
                    else if((pbound_Y==1 || resbound_Y==1) && j==sizeY2-1-isoNum && threadj==numPy-1){
                        countPyOut=countPyOut+1;
                        listPyOutnodes[countPyOut]=tempID;
                    }
                    else if((pbound_Z==1 || resbound_Z==1) && k==isoNum && threadk==0){
                        countPzIn=countPzIn+1;
                        listPzInnodes[countPzIn]=tempID;
                    }
                    else if((pbound_Z==1 || resbound_Z==1) && k==sizeZ2-1-isoNum && threadk==numPz-1){
                        countPzOut=countPzOut+1;
                        listPzOutnodes[countPzOut]=tempID;
                    }
                    else{
                        countFint=countFint+1;
                        listFintnodes[countFint]=tempID;
                    }
                }
                if (bbRegionSub[tempID]==1){
                    if((pbound_X==1 || resbound_X==1) && i==isoNum && threadi==0){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else if((pbound_X==1 || resbound_X==1) && i==sizeX2-1-isoNum && threadi==numPx-1){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else if((pbound_Y==1 || resbound_Y==1) && j==isoNum && threadj==0){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else if((pbound_Y==1 || resbound_Y==1) && j==sizeY2-1-isoNum && threadj==numPy-1){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else if((pbound_Z==1 || resbound_Z==1) && k==isoNum && threadk==0){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else if((pbound_Z==1 || resbound_Z==1) && k==sizeZ2-1-isoNum && threadk==numPz-1){
                        countBBP=countBBP+1;
                        listBBPnodes[countBBP]=tempID;
                    }
                    else{
                        countBB = countBB+1;
                        listBBnodes[countBB]=tempID;
                    }
                }
            }
        }
    }

    ///output number of bounce-back and fluid nodes
    int countBB_sum=0;
    int countF_sum=0;
    int countBBP_sum=0;
    int countBBsend, countFsend;
    if (rank>0){
        countBBsend=countBB+1;
        countFsend=countF+1;
        MPI_Send(&countBBsend, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);  //send all fluid node counts to rank 0, tag=2
        MPI_Send(&countFsend, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);  //send all bb node counts to rank 0, tag=3
    }
    else{  //rank 0 collects all other processor fluid and bb node counts
        countBB_sum=countBB+1;
        countF_sum=countF+1;
        for (thread=1;thread<=size-1;thread++){
            MPI_Recv(&countBBsend, 1, MPI_INT, thread, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            countBB_sum=countBB_sum+countBBsend;
            MPI_Recv(&countFsend, 1, MPI_INT, thread, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            countF_sum=countF_sum+countFsend;
        }
        printf("BB nodes (not including BB nodes on pressure/reservoir boundaries) in lattice:  %d\nFluid nodes in lattice: %d\n", countBB_sum, countF_sum);
    }

///+++++++++++++++Print out local mean or global viscosity model+++++++++++++++++////
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///

    ///Set to global viscosity model for standard simulation
    //collect individual sv and sq
    int svCount_sum=0;
    int sqCount_sum=0;
    double svSum_sum=0;
    double sqSum_sum=0;

    if(globalVisc_==1){
        mu=nu*rhoFluid;
        if(rank==0)printf("Global Viscosity Model (Standard) mean kinematic visc = %f, mean dynamic visc = %f\n", nu,mu);
    }
    else{
        if (rank>0){
            MPI_Send(&svCount, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
            MPI_Send(&sqCount, 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
            MPI_Send(&svSum, 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
            MPI_Send(&sqSum, 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
        }
        else{  //rank 0 collects all other processor sv and sq
            svCount_sum=svCount;
            sqCount_sum=sqCount;
            svSum_sum=svSum;
            sqSum_sum=sqSum;
            for (thread=1;thread<=size-1;thread++){
                MPI_Recv(&svCount, 1, MPI_DOUBLE, thread, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                svCount_sum=svCount_sum+svCount;
                MPI_Recv(&sqCount, 1, MPI_DOUBLE, thread, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sqCount_sum=sqCount_sum+sqCount;
                MPI_Recv(&svSum, 1, MPI_DOUBLE, thread, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                svSum_sum=svSum_sum+svSum;
                MPI_Recv(&sqSum, 1, MPI_DOUBLE, thread, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sqSum_sum=sqSum_sum+sqSum;
            }
            svMean=svSum_sum/svCount_sum;
            sqMean=sqSum_sum/sqCount_sum;
            svUnbound=1./(0.5+pow((6./pi),0.5)*mfp);
            nu = (cs*cs)*((1./svUnbound)-0.5);
            mu = nu*rhoFluid;
        }
        if(rank==0){
            if(stdBB_==0){printf("Second-order N-S Slip Flow with CDBB-Landry Model, mean sv = %f\n mean sq = %f\n mean unbound kinematic visc = %f\n mean unbound dynamic visc = %f\n", svMean, sqMean, nu, mu);}
            else{printf("StandardBB-LEV Model, why are you using this setup? is this a mistake? mean sv = %f\n mean sq = %f\n mean unbound kinematic visc = %f\n mean unbound dynamic visc = %f\n", svMean, sqMean, nu, mu);}
        }
    }
///+++++++++++++++Begin initialization of lattice particle distributions and macroscopic vectors+++++++++///
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///
    ///Initialize the fluid nodes with random density, this can be superseded by pressure gradient or read in density file
    //fill in rho with rhoWall on solid nodes, vel and population arrays with zeros
    for(k=0;k<sizeZ2;k++){
        for(j=0;j<sizeY2;j++){
            for(i=0;i<sizeX2;i++){
                tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                if(bbRegionSub[tempID]==0) {	//0->pore, 1->boundary
                    rho[tempID]=rhoFluid+rhoFluidran*(-1.+2*(double)(rand()%1000)/1000.);
                    for(pop=0;pop<=18;pop++){
                        tempIDpop=(sizepop*tempID)+pop;
                        fOut[tempIDpop]=t[pop]*rho[tempID];
                    }
                }
                else {rho[tempID]=rhoWall;}
                ux[tempID]=0.;
                uy[tempID]=0.;
                uz[tempID]=0.;
                forceAllx[tempID]=0.;
                forceAlly[tempID]=0.;
                forceAllz[tempID]=0.;
                for(pop=0;pop<19;pop++){
                    tempIDpop=(sizepop*tempID)+pop;
                    fIn[tempIDpop]=0.;
                }
            }
        }
    }

///++++++++++++Linear pressure distribution intialization+++++++++++++++++++++++++++++////
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++////

//pressure gradients and/or forceExts to be used in permeability determination
//If pressure boundaries exist then densities and initial distribution function should reflect pressure gradient, if no pressure gradient then initialize with rhoFluid
    //Determine subdomains
    int subx = floor(lx/numPx);
    int suby = floor(ly/numPy);
    int subz = floor(lz/numPz);

    int a;
    int xstart[numPx];
    int xend[numPx];
    xstart[0]=0;
    xend[numPx-1]=lx-1;
    for(a=1;a<=numPx-1;a++){
        xstart[a]=subx*a;
        xend[a-1]=(subx*a)-1;
    }
    int ystart[numPy];
    int yend[numPy];
    ystart[0]=0;
    yend[numPy-1]=ly-1;
    for(a=1;a<=numPy-1;a++){
        ystart[a]=suby*a;
        yend[a-1]=(suby*a)-1;
    }
    int zstart[numPz];
    int zend[numPz];
    zstart[0]=0;
    zend[numPz-1]=lz-1;
    for(a=1;a<=numPz-1;a++){
        zstart[a]=subz*a;
        zend[a-1]=(subz*a)-1;
    }

    double gradPX, gradPY, gradPZ;
    if (periodic_x==1 && accell[0]>0.0){
        gradPX = rhoFluid*(double)lx*(accell[0])/(double)lx;
        flowDir=0;
        gradP_perm=gradPX;
        if(rank==0)printf("Force Pressure Gradient in x: %E\n", gradPX);
    }
    if (periodic_y==1 && accell[1]>0.0){
        gradPY = rhoFluid*(double)ly*(accell[1])/(double)ly;
        flowDir=1;
        gradP_perm=gradPY;
        if(rank==0)printf("Force Pressure Gradient in y: %E\n", gradPY);
    }
    if (periodic_z==1 && accell[2]>0.0){
        gradPZ = rhoFluid*(double)lz*(accell[2])/(double)lz;
        flowDir=2;
        gradP_perm=gradPZ;
        if(rank==0)printf("Force Pressure Gradient in z: %E\n", gradPZ);
    }

    if (pbound_X==1){
        gradPX=(pInX+pOutX)/((double)lx-1.0);
        gradP_perm=gradPX;
        flowDir=0;
        rhoFluid=3.*((pInX+pOutX)/2);
        mu=nu*rhoFluid;
        if(rank==0){
            printf("Pressure Gradient in x: %E\n", gradPX);
        }
        if (pbound_X_LinP==1){
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp = (pInX*3.) - (xstart[threadi]+i-1)*(((pInX*3.)-(pOutX*3.))/((double)lx-1));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fOut[tempIDpop]=rhoTemp*t[pop];
                }
            }
        }
    }
    if (pbound_Y==1){
        gradPY=(pInY-pOutY)/((double)ly-1.0);
        gradP_perm=gradPY;
        flowDir=1;
        rhoFluid=3.*((pInY+pOutY)/2);
        mu=nu*rhoFluid;
        if(rank==0){
            printf("Pressure Gradient in y: %E\n", gradPY);
        }
        if (pbound_Y_LinP==1){
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp = (pInY*3.) - (ystart[threadj]+j-1)*(((pInY*3.)-(pOutY*3.))/((double)ly-1.0));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fOut[tempIDpop]=rhoTemp*t[pop];
                }
            }
        }
    }
    if (pbound_Z==1){
        gradPZ=(pInZ-pOutZ)/((double)lz-1.0);
        gradP_perm=gradPZ;
        flowDir=2;
        rhoFluid=3.*((pInZ+pOutZ)/2);
        mu=nu*rhoFluid;
        if(rank==0){
            printf("Pressure Gradient in z: %E\n", gradPZ);
        }
        if (pbound_Z_LinP==1){
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp = (pInZ*3.) - (zstart[threadk]+k-1)*(((pInZ*3.)-(pOutZ*3.))/((double)lz-1.0));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fOut[tempIDpop]=rhoTemp*t[pop];
                }
            }
        }
    }
    if(rank==0) printf("Pressure Gradient used in permeability calculation: %E\n", gradP_perm);

///+++++++++++++++++++++++ Read in rho file initialization +++++++++++++++++++++++++++++++///
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++///


    ///Pass indexFunc between processes, indexfunction indicates node type 1 for solid and BB, 0 for fluid
    for(kk=0;kk<numCommR;kk++){
        eDir=passDirR[kk];
        countS=-1;
        for(k=zRSendBounds[eDir][0];k<=zRSendBounds[eDir][1];k++){
            for(j=yRSendBounds[eDir][0];j<=yRSendBounds[eDir][1];j++){
                for(i=xRSendBounds[eDir][0];i<=xRSendBounds[eDir][1];i++){
                    tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                    ind_send[++countS]=indexFunc[tempID];
                }
            }
        }
        tagS = (3*size+rank+1)*100 + eDir; //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
        commRType=oppComm[eDir];
        tagR = (3*size+commRank[eDir]+1)*100 + commRType; //  //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
        MPI_Sendrecv(ind_send, sizeRSend[eDir], MPI_INT, commRank[eDir], tagS, ind_recv, sizeRSend[eDir], MPI_DOUBLE, commRank[eDir], tagR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        countS=-1;
        for(k=zRRecvBounds[commRType][0];k<=zRRecvBounds[commRType][1];k++){
            for(j=yRRecvBounds[commRType][0];j<=yRRecvBounds[commRType][1];j++){
                for(i=xRRecvBounds[commRType][0];i<=xRRecvBounds[commRType][1];i++){
                    tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                    indexFunc[tempID]=ind_recv[++countS];
                }
            }
        }
    }

    ///read in rho, if rho is read in it will supersede default random density or option linear gradient
    char rhoFileO[350];
    if(readRhoIn==1){
        sprintf(rhoFileO, "%s/%s_%d", rhoFile, rhoFile, rank);
        FILE *rhoIn = fopen(rhoFileO,"r");
        if (rhoIn == NULL) {
            perror("error opening density file");
        }
        for(k=isoNum;k<sizeZ2-isoNum;k++){
            for(j=isoNum;j<sizeY2-isoNum;j++){
                for(i=isoNum;i<sizeX2-isoNum;i++){
                    tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                    fscanf(rhoIn,"%lf",&rho[tempID]);
                    rho_sum=rho_sum+rho[tempID];
                    for(pop=0;pop<=18;pop++){
                        tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                        fOut[tempIDpop]=rho[tempID]*t[pop]*(1.-(double)indexFunc[tempID]); //only fill in fOut for fluid nodes
                    }
                }
            }
        }
        fclose(rhoIn);
    }

///++++++++++++++++++++++++++++++++++ MAIN LOOP ++++++++++++++++++++++++++++++++++++++++++++//
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    iT=0;
    while(iT<maxT){
        iT=iT+1;
        previT=iT;

        //STREAMING STEP, all fluid nodes are streamed
        for(iD=0; iD<=countF; iD++){
            tempID=listFnodes[iD];
            k=floor(tempID/((sizeX2)*(sizeY2)));
            j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
            i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
            for(pop=0;pop<=18;pop++){
                tempIDpop=(sizepop*tempID)+pop;
                tempIDpop2=((k+cz[pop])*sizepop*sizeX2*sizeY2)+((j+cy[pop])*sizepop*sizeX2)+((i+cx[pop])*sizepop)+pop;
                fIn[tempIDpop2]=fOut[tempIDpop];
            }
        }


///+++++++++++++++++++++++++++++++++++++++++++++STREAMING MPI COMMUNICATION+++++++++++++++++++++++++++++++++++//
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        for(kk=0;kk<numComm;kk++){
            eDir=passDir[kk];
            countS=-1;
            for(k=zSendBounds[eDir][0];k<=zSendBounds[eDir][1];k++){
                for(j=ySendBounds[eDir][0];j<=ySendBounds[eDir][1];j++){
                    for(i=xSendBounds[eDir][0];i<=xSendBounds[eDir][1];i++){
                        for(pop=0;pop<numPop[eDir];pop++){
                            popRS=passPops[eDir][pop];
                            tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+popRS;
                            fIn_send[++countS]=fIn[tempIDpop];
                        }
                    }
                }
            }
            tagS = (rank+1)*100 + eDir; //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            commRType=oppComm[eDir];
            tagR = (commRank[eDir]+1)*100 + commRType; //  //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            //printf("rank = %i, eDir = %i, tagS = %i, tagR = %i, countS = %i, sizeSend = %i\n", rank,eDir,tagS,tagR,countS,sizeSend[eDir]);
            MPI_Sendrecv(fIn_send, sizeSend[eDir], MPI_DOUBLE, commRank[eDir], tagS, fIn_recv, sizeSend[eDir], MPI_DOUBLE, commRank[eDir], tagR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            countS=-1;
            for(k=zRecvBounds[commRType][0];k<=zRecvBounds[commRType][1];k++){
                for(j=yRecvBounds[commRType][0];j<=yRecvBounds[commRType][1];j++){
                    for(i=xRecvBounds[commRType][0];i<=xRecvBounds[commRType][1];i++){
                        for(pop=0;pop<numPop[eDir];pop++){
                            popRS=passPops[commRType][pop];
                            tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+popRS;
                            fIn[tempIDpop]=fIn_recv[++countS];
                        }
                    }
                }
            }
        }



///+++++++++++++++++++++++++++++++++++++ Bounce-Back Schemes ++++++++++++++++++++++++++++++++++++++++++//
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        ///standard halfway bb
        if(stdBB_==1){
            for(iD=0; iD<=countBB; iD++){
                tempID=listBBnodes[iD];		// list BB nodes
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=1;pop<=18;pop++){
                    popRS=opp[pop];
                    tempIDpop=(sizepop*tempID)+popRS;
                    tempIDpop2=((k+cz[pop])*sizepop*sizeX2*sizeY2)+((j+cy[pop])*sizepop*sizeX2)+((i+cx[pop])*sizepop)+pop;
                    fIn[tempIDpop2]=fIn[tempIDpop];
                }
            }
            //BounceBack step on pressure boundary BB nodes is only standard halfway BB
            for(iD=0; iD<=countBBP; iD++){
                tempID=listBBPnodes[iD];	//list BBP nodes, pressure boundary
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=1;pop<=18;pop++){
                    popRS=opp[pop];
                    tempIDpop=(sizepop*tempID)+popRS;
                    tempIDpop2=((k+cz[pop])*sizepop*sizeX2*sizeY2)+((j+cy[pop])*sizepop*sizeX2)+((i+cx[pop])*sizepop)+pop;
                    fIn[tempIDpop2]=fIn[tempIDpop];
                }
            }
        }
        else {
            ///combined diffusive standard bounceback
            for(iD=0; iD<=countBB; iD++){
                tempID=listBBnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                An=0;
                fSum=0;
                for(pop=1;pop<=18;pop++){		//
                    tempIDpop=(sizepop*tempID)+pop;
                    An=An+(t[pop]*(double)fluidCommSub[tempIDpop]);	//sum of the D3Q19 weights (wi) of the incoming fi.
                    fSum=fSum+fIn[tempIDpop]*(double)fluidCommSub[tempIDpop];	//sum of the incoming fi for the solid boundary node
                }
                sigma=fSum/An;
                for(pop=1;pop<=18;pop++){
                    popRS=opp[pop];
                    tempIDpop=(sizepop*tempID)+popRS;
                    tempIDpop2=((k+cz[pop])*sizepop*sizeX2*sizeY2)+((j+cy[pop])*sizepop*sizeX2)+((i+cx[pop])*sizepop)+pop;
                    fIn[tempIDpop2]=
                        ((rBB*(sigma*t[pop]*uWallEq[pop])+(1.-rBB)*fIn[tempIDpop])*(double)fluidCommSub[tempIDpop])
                        +(fIn[tempIDpop2]*(1.-(double)fluidCommSub[tempIDpop])); //without this line fIn[tempIDpop2] is set equal to 0
                }
            }
            //BounceBack step on pressure boundary BB nodes is only standard halfway BB
            for(iD=0; iD<=countBBP; iD++){
                tempID=listBBPnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=1;pop<=18;pop++){
                    popRS=opp[pop];
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+popRS;
                    tempIDpop2=((k+cz[pop])*sizepop*sizeX2*sizeY2)+((j+cy[pop])*sizepop*sizeX2)+((i+cx[pop])*sizepop)+pop;
                    fIn[tempIDpop2]=fIn[tempIDpop];
                }
            }
        }



///+++++++++++++++++++++++++++++++++++++ START BOUNCEBACK MPI COMMUNICATION ++++++++++++++++++++++++++++++++++++++//
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        for(kk=0;kk<numComm;kk++){
            eDir=passDir[kk];
            countS=-1;
            for(k=zSendBounds[eDir][0];k<=zSendBounds[eDir][1];k++){
                for(j=ySendBounds[eDir][0];j<=ySendBounds[eDir][1];j++){
                    for(i=xSendBounds[eDir][0];i<=xSendBounds[eDir][1];i++){
                        for(pop=0;pop<numPop[eDir];pop++){
                            popRS=passPops[eDir][pop];
                            tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+popRS;
                            fIn_send[++countS]=fIn[tempIDpop];
                        }
                    }
                }
            }
            tagS = (size+rank+1)*100 + eDir; //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            commRType=oppComm[eDir];
            tagR = (size+commRank[eDir]+1)*100 + commRType; //  //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            MPI_Sendrecv(fIn_send, sizeSend[eDir], MPI_DOUBLE, commRank[eDir], tagS, fIn_recv, sizeSend[eDir], MPI_DOUBLE, commRank[eDir], tagR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            countS=-1;
            for(k=zRecvBounds[commRType][0];k<=zRecvBounds[commRType][1];k++){
                for(j=yRecvBounds[commRType][0];j<=yRecvBounds[commRType][1];j++){
                    for(i=xRecvBounds[commRType][0];i<=xRecvBounds[commRType][1];i++){
                        for(pop=0;pop<numPop[eDir];pop++){
                            popRS=passPops[commRType][pop];
                            tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+popRS;
                            fIn[tempIDpop]=fIn_recv[++countS];
                        }
                    }
                }
            }
        }

///+++++++++++++++++++++++++++++++++++++ END BOUNCEBACK MPI COMMUNICATION ++++++++++++++++++++++++++++++++++++++++//
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        //If pressure boundaries exist
        //Using Hecht Harting 2010 boundary conditions

        if(pbound_Z==1 && threadk==0){
            for(iD=0; iD<=countPzIn; iD++){
                tempID=listPzInnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                uz[tempID]=1.-(1./(pInZ*3.))*(fCalc[0]+fCalc[1]+fCalc[2]+fCalc[3]+fCalc[4]+fCalc[7]+fCalc[9]+fCalc[8]+fCalc[10]+2.*(fCalc[6]+fCalc[13]+fCalc[14]+fCalc[17]+fCalc[18]));
                NxZ=(1./2.)*(fCalc[1]+fCalc[7]+fCalc[9]-(fCalc[2]+fCalc[8]+fCalc[10]))-(1./3.)*(pInZ*3.)*pInZ_ux;
                NyZ=(1./2.)*(fCalc[3]+fCalc[7]+fCalc[8]-(fCalc[4]+fCalc[9]+fCalc[10]))-(1./3.)*(pInZ*3.)*pInZ_uy;
                fCalc[5]=fCalc[6]+(1./3.)*(pInZ*3.)*uz[tempID];
                fCalc[11]=fCalc[14]+(1./6.)*(pInZ*3.)*(uz[tempID]+pInZ_ux)-NxZ;
                fCalc[12]=fCalc[13]+(1./6.)*(pInZ*3.)*(uz[tempID]-pInZ_ux)+NxZ;
                fCalc[15]=fCalc[18]+(1./6.)*(pInZ*3.)*(uz[tempID]+pInZ_uy)-NyZ;
                fCalc[16]=fCalc[17]+(1./6.)*(pInZ*3.)*(uz[tempID]-pInZ_uy)+NyZ;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+5;
                fIn[tempIDpop]=fCalc[5];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+11;
                fIn[tempIDpop]=fCalc[11];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+12;
                fIn[tempIDpop]=fCalc[12];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+15;
                fIn[tempIDpop]=fCalc[15];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+16;
                fIn[tempIDpop]=fCalc[16];
            }
        }
        if(pbound_Z==1 && threadk==numPz-1){
            for(iD=0; iD<=countPzOut; iD++){
                tempID=listPzOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                uz[tempID]=-1.+(1./(pOutZ*3.))*(fCalc[0]+fCalc[1]+fCalc[2]+fCalc[3]+fCalc[4]+fCalc[7]+fCalc[9]+fCalc[8]+fCalc[10]+2.*(fCalc[5]+fCalc[11]+fCalc[12]+fCalc[15]+fCalc[16]));
                NxZ=(1./2.)*(fCalc[1]+fCalc[7]+fCalc[9]-(fCalc[2]+fCalc[8]+fCalc[10]))-(1./3.)*(pOutZ*3.)*pOutZ_ux;
                NyZ=(1./2.)*(fCalc[3]+fCalc[7]+fCalc[8]-(fCalc[4]+fCalc[9]+fCalc[10]))-(1./3.)*(pOutZ*3.)*pOutZ_uy;
                fCalc[6]=fCalc[5]-(1./3.)*(pOutZ*3.)*uz[tempID];
                fCalc[13]=fCalc[12]+(1./6.)*(pOutZ*3.)*(-uz[tempID]+pOutZ_ux)-NxZ;
                fCalc[14]=fCalc[11]+(1./6.)*(pOutZ*3.)*(-uz[tempID]-pOutZ_ux)+NxZ;
                fCalc[17]=fCalc[16]+(1./6.)*(pOutZ*3.)*(-uz[tempID]+pOutZ_uy)-NyZ;
                fCalc[18]=fCalc[15]+(1./6.)*(pOutZ*3.)*(-uz[tempID]-pOutZ_uy)+NyZ;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+6;
                fIn[tempIDpop]=fCalc[6];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+13;
                fIn[tempIDpop]=fCalc[13];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+14;
                fIn[tempIDpop]=fCalc[14];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+17;
                fIn[tempIDpop]=fCalc[17];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+18;
                fIn[tempIDpop]=fCalc[18];
            }
        }

        if(pbound_Y==1 && threadj==0){
            for(iD=0; iD<=countPyIn; iD++){
                tempID=listPyInnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                uy[tempID]=1.-(1./(pInY*3.))*(fCalc[0]+fCalc[1]+fCalc[2]+fCalc[5]+fCalc[6]+fCalc[11]+fCalc[13]+fCalc[12]+fCalc[14]+2.*(fCalc[4]+fCalc[9]+fCalc[10]+fCalc[16]+fCalc[18]));
                NxY=(1./2.)*(fCalc[1]+fCalc[11]+fCalc[13]-(fCalc[2]+fCalc[12]+fCalc[14]))-(1./3.)*(pInY*3.)*pInY_ux;
                NzY=(1./2.)*(fCalc[5]+fCalc[11]+fCalc[12]-(fCalc[6]+fCalc[13]+fCalc[14]))-(1./3.)*(pInY*3.)*pInY_uz;
                fCalc[3]=fCalc[4]+(1./3.)*(pInY*3.)*uy[tempID];
                fCalc[7]=fCalc[10]+(1./6.)*(pInY*3.)*(uy[tempID]+pInY_ux)-NxY;
                fCalc[8]=fCalc[9]+(1./6.)*(pInY*3.)*(uy[tempID]-pInY_ux)+NxY;
                fCalc[15]=fCalc[18]+(1./6.)*(pInY*3.)*(uy[tempID]+pInY_uz)-NzY;
                fCalc[17]=fCalc[16]+(1./6.)*(pInY*3.)*(uy[tempID]-pInY_uz)+NzY;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+3;
                fIn[tempIDpop]=fCalc[3];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+7;
                fIn[tempIDpop]=fCalc[7];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+8;
                fIn[tempIDpop]=fCalc[8];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+15;
                fIn[tempIDpop]=fCalc[15];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+17;
                fIn[tempIDpop]=fCalc[17];
            }
        }
        if(pbound_Y==1 && threadj==numPy-1){
            for(iD=0; iD<=countPyOut; iD++){
                tempID=listPyOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                uy[tempID]=-1.+(1./(pOutY*3.))*(fCalc[0]+fCalc[1]+fCalc[2]+fCalc[5]+fCalc[6]+fCalc[11]+fCalc[13]+fCalc[12]+fCalc[14]+2.*(fCalc[3]+fCalc[7]+fCalc[8]+fCalc[15]+fCalc[17]));
                NxY=(1./2.)*(fCalc[1]+fCalc[11]+fCalc[13]-(fCalc[2]+fCalc[12]+fCalc[14]))-(1./3.)*(pInY*3.)*pOutY_ux;
                NzY=(1./2.)*(fCalc[5]+fCalc[11]+fCalc[12]-(fCalc[6]+fCalc[13]+fCalc[14]))-(1./3.)*(pInY*3.)*pOutY_uz;
                fCalc[4]=fCalc[3]-(1./3.)*(pOutY*3.)*uy[tempID];
                fCalc[10]=fCalc[7]+(1./6.)*(pOutY*3.)*(-uy[tempID]-pOutY_ux)+NxY;
                fCalc[9]=fCalc[8]+(1./6.)*(pOutY*3.)*(-uy[tempID]+pOutY_ux)-NxY;
                fCalc[18]=fCalc[15]+(1./6.)*(pOutY*3.)*(-uy[tempID]-pOutY_uz)+NzY;
                fCalc[16]=fCalc[17]+(1./6.)*(pOutY*3.)*(-uy[tempID]+pOutY_uz)-NzY;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+4;
                fIn[tempIDpop]=fCalc[4];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+10;
                fIn[tempIDpop]=fCalc[10];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+9;
                fIn[tempIDpop]=fCalc[9];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+18;
                fIn[tempIDpop]=fCalc[18];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+16;
                fIn[tempIDpop]=fCalc[16];
            }
        }

        if(pbound_X==1 && threadi==0){
            for(iD=0; iD<=countPxIn; iD++){
                tempID=listPxInnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                ux[tempID]=1.-(1./(pInX*3.))*(fCalc[0]+fCalc[3]+fCalc[4]+fCalc[5]+fCalc[6]+fCalc[15]+fCalc[17]+fCalc[16]+fCalc[18]+2.*(fCalc[2]+fCalc[8]+fCalc[10]+fCalc[12]+fCalc[14]));
                NyX=(1./2.)*(fCalc[3]+fCalc[15]+fCalc[17]-(fCalc[4]+fCalc[16]+fCalc[18]))-(1./3.)*(pInX*3.)*pInX_uy;
                NzX=(1./2.)*(fCalc[5]+fCalc[15]+fCalc[16]-(fCalc[6]+fCalc[17]+fCalc[18]))-(1./3.)*(pInX*3.)*pInX_uz;
                fCalc[1]=fCalc[2]+(1./3.)*(pInX*3.)*ux[tempID];
                fCalc[7]=fCalc[10]+(1./6.)*(pInX*3.)*(ux[tempID]+pInX_uy)-NyX;
                fCalc[9]=fCalc[8]+(1./6.)*(pInX*3.)*(ux[tempID]-pInX_uy)+NyX;
                fCalc[11]=fCalc[14]+(1./6.)*(pInX*3.)*(ux[tempID]+pInX_uz)-NzX;
                fCalc[13]=fCalc[12]+(1./6.)*(pInX*3.)*(ux[tempID]-pInX_uz)+NzX;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+1;
                fIn[tempIDpop]=fCalc[1];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+7;
                fIn[tempIDpop]=fCalc[7];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+9;
                fIn[tempIDpop]=fCalc[9];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+11;
                fIn[tempIDpop]=fCalc[11];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+13;
                fIn[tempIDpop]=fCalc[13];
            }
        }
        if(pbound_X==1 && threadi==numPx-1){
            for(iD=0; iD<=countPxOut; iD++){
                tempID=listPxOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0;pop<=18;pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fCalc[pop]=fIn[tempIDpop];
                }
                ux[tempID]=-1.+(1./(pOutX*3.))*(fCalc[0]+fCalc[3]+fCalc[4]+fCalc[5]+fCalc[6]+fCalc[15]+fCalc[17]+fCalc[16]+fCalc[18]+2.*(fCalc[1]+fCalc[7]+fCalc[9]+fCalc[11]+fCalc[13]));
                NyX=(1./2.)*(fCalc[3]+fCalc[15]+fCalc[17]-(fCalc[4]+fCalc[16]+fCalc[18]))-(1./3.)*(pInX*3.)*pOutX_uy;
                NzX=(1./2.)*(fCalc[5]+fCalc[15]+fCalc[16]-(fCalc[6]+fCalc[17]+fCalc[18]))-(1./3.)*(pInX*3.)*pOutX_uz;
                fCalc[2]=fCalc[1]-(1./3.)*(pOutX*3.)*ux[tempID];
                fCalc[8]=fCalc[9]+(1./6.)*(pOutX*3.)*(-ux[tempID]+pOutX_uy)-NyX;
                fCalc[10]=fCalc[7]+(1./6.)*(pOutX*3.)*(-ux[tempID]-pOutX_uy)+NyX;
                fCalc[12]=fCalc[13]+(1./6.)*(pOutX*3.)*(-ux[tempID]+pOutX_uz)-NzX;
                fCalc[14]=fCalc[11]+(1./6.)*(pOutX*3.)*(-ux[tempID]-pOutX_uz)+NzX;
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+2;
                fIn[tempIDpop]=fCalc[2];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+8;
                fIn[tempIDpop]=fCalc[8];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+10;
                fIn[tempIDpop]=fCalc[10];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+12;
                fIn[tempIDpop]=fCalc[12];
                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+14;
                fIn[tempIDpop]=fCalc[14];
            }
        }

        if(resbound_Z==1 && threadk==0){
            for(iD=0; iD<=countPzIn; iD++){
                tempID=listPzInnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoInZ;
                }
            }
        }
        if(resbound_Z==1 && threadk==numPz-1){
            for(iD=0; iD<=countPzOut; iD++){
                tempID=listPzOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoOutZ;
                }
            }
        }
        if(resbound_Y==1 && threadj==0){
            for(iD=0; iD<=countPyIn; iD++){
                tempID=listPyInnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoInY;
                }
            }
        }
        if(resbound_Y==1 && threadj==numPz-1){
            for(iD=0; iD<=countPyOut; iD++){
                tempID=listPyOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoOutY;
                }
            }
        }
        if(resbound_X==1 && threadi==0){
            for(iD=0; iD<=countPxIn; iD++){
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoInX;
                }
            }
        }
        if(resbound_X==1 && threadi==numPz-1){
            for(iD=0; iD<=countPxOut; iD++){
                tempID=listPxOutnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                for(pop=0; pop<=18; pop++){
                    tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                    fIn[tempIDpop]=t[pop]*resRhoOutX;
                }
            }
        }


///++++++++++++++++++++++++++++++++ Calculate Macroscopic variables and Collide ++++++++++++++++++++++++++++//
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        ///Rho calculated on all fluid nodes
        for(iD=0; iD<=countF; iD++){
            tempID=listFnodes[iD];
            rhoTemp=0.;
            for(pop=0;pop<=18;pop++){
                tempIDpop=sizepop*tempID+pop;
                rhoTemp=rhoTemp+fIn[tempIDpop];
            }
            ///unroll for loop option, suppose to be faster? maybe not.
            //rho[tempID]==fIn[tempIDpop+0]+fIn[tempIDpop+1]+fIn[tempIDpop+2]+fIn[tempIDpop+3]+fIn[tempIDpop+4]+
            //        fIn[tempIDpop+5]+fIn[tempIDpop+6]+fIn[tempIDpop+7]+fIn[tempIDpop+8]+fIn[tempIDpop+9]+
            //        fIn[tempIDpop+10]+fIn[tempIDpop+11]+fIn[tempIDpop+12]+fIn[tempIDpop+13]+fIn[tempIDpop+14]+
            //        fIn[tempIDpop+15]+fIn[tempIDpop+16]+fIn[tempIDpop+17]+fIn[tempIDpop+18];
            rho[tempID]=rhoTemp;
        }
        ///pass Rho to ghost nodes SC forcing
        for(kk=0;kk<numCommR;kk++){
            eDir=passDirR[kk];
            countS=-1;
            for(k=zRSendBounds[eDir][0];k<=zRSendBounds[eDir][1];k++){
                for(j=yRSendBounds[eDir][0];j<=yRSendBounds[eDir][1];j++){
                    for(i=xRSendBounds[eDir][0];i<=xRSendBounds[eDir][1];i++){
                        tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                        scalar_send[++countS]=rho[tempID];
                    }
                }
            }
            tagS = (2*size+rank+1)*100 + eDir; //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            commRType=oppComm[eDir];
            tagR = (2*size+commRank[eDir]+1)*100 + commRType; //  //tag id coding, as in tag = 101 means from thread 1 in the direction of e1 set equal to 1
            MPI_Sendrecv(scalar_send, sizeRSend[eDir], MPI_DOUBLE, commRank[eDir], tagS, scalar_recv, sizeRSend[eDir], MPI_DOUBLE, commRank[eDir], tagR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            countS=-1;
            for(k=zRRecvBounds[commRType][0];k<=zRRecvBounds[commRType][1];k++){
                for(j=yRRecvBounds[commRType][0];j<=yRRecvBounds[commRType][1];j++){
                    for(i=xRRecvBounds[commRType][0];i<=xRRecvBounds[commRType][1];i++){
                        tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                        rho[tempID]=scalar_recv[++countS];
                    }
                }
            }
        }

        ///Shan-Chen SCMP pseudo-potential model, not applied to pressure/reservoir boundaries
        if(isoNum==1){
            for(iD=0; iD<=countFint; iD++){
                tempID=listFintnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp1=rho[tempID];
                psiSC1=rhoSC*(1.-exp(0.-rhoTemp1/rhoSC)); // ShanChen Original EOS, 1993
                forceSCx=0.;
                forceSCy=0.;
                forceSCz=0.;
                for(pop=1;pop<=18;pop++){
                    tempID2=((k+cz[pop])*sizeX2*sizeY2)+((j+cy[pop])*sizeX2)+(i+cx[pop]);
                    rhoTemp2=rho[tempID2];
                    psiSC2=rhoSC*(1.-exp(0.-rhoTemp2/rhoSC)); //fluid-solid interaction strength determined by rhowall
                    forceSCx=forceSCx-(cs*cs*potGSC*weightSC[pop]*psiSC1*psiSC2*cxD[pop]);
                    forceSCy=forceSCy-(cs*cs*potGSC*weightSC[pop]*psiSC1*psiSC2*cyD[pop]);
                    forceSCz=forceSCz-(cs*cs*potGSC*weightSC[pop]*psiSC1*psiSC2*czD[pop]);
                }
                //sum interfacial and exterior forces
                forceAllx[tempID]=forceSCx+accell[0]*rhoTemp1;
                forceAlly[tempID]=forceSCy+accell[1]*rhoTemp1;
                forceAllz[tempID]=forceSCz+accell[2]*rhoTemp1;
                //if(rank==14 && i==50 && k==1)printf("iT=%d,i=%d,j=%d,k=%d forces: %E %E %E psi1 = %E, psi2e18 = %E\n",iT,i,j,k,forceAll[0][i][j][k],forceAll[1][i][j][k],forceAll[2][i][j][k], psiSC1, psiSC2);
            }
        }
        else if (isoNum==2){
            for(iD=0; iD<=countFint; iD++){
                tempID=listFintnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp1=rho[tempID];
                psiSC1=rhoSC*(1.-exp(0.-rhoTemp1/rhoSC)); // ShanChen Original EOS
                forceSCx=0.;
                forceSCy=0.;
                forceSCz=0.;
                for(pop=0;pop<92;pop++){
                    tempID2=((k+cz8[pop])*sizeX2*sizeY2)+((j+cy8[pop])*sizeX2)+(i+cx8[pop]);
                    rhoTemp2=rho[tempID2];
                    psiSC2=rhoSC*(1.-exp(0.-rhoTemp2/rhoSC)); //fluid-solid interaction strength determined by rhowall
                    forceSCx=forceSCx-(cs*cs*potGSC*weight8[pop]*psiSC1*psiSC2*(double)cx8[pop]);
                    forceSCy=forceSCy-(cs*cs*potGSC*weight8[pop]*psiSC1*psiSC2*(double)cy8[pop]);
                    forceSCz=forceSCz-(cs*cs*potGSC*weight8[pop]*psiSC1*psiSC2*(double)cz8[pop]);
                }
                //sum interfacial and exterior forces
                forceAllx[tempID]=forceSCx+accell[0]*rhoTemp1;
                forceAlly[tempID]=forceSCy+accell[1]*rhoTemp1;
                forceAllz[tempID]=forceSCz+accell[2]*rhoTemp1;
                //if(rank==14 && i==50 && k==1)printf("iT=%d,i=%d,j=%d,k=%d forces: %E %E %E psi1 = %E, psi2e18 = %E\n",iT,i,j,k,forceAll[0][i][j][k],forceAll[1][i][j][k],forceAll[2][i][j][k], psiSC1, psiSC2);
            }
        }
        else{
            for(iD=0; iD<=countFint; iD++){
                tempID=listFintnodes[iD];
                k=floor(tempID/((sizeX2)*(sizeY2)));
                j=floor((tempID-(k*(sizeX2)*(sizeY2)))/(sizeX2));
                i=tempID-(j*(sizeX2))-(k*(sizeX2)*(sizeY2));
                rhoTemp1=rho[tempID];
                psiSC1=rhoSC*(1.-exp(0.-rhoTemp1/rhoSC)); // ShanChen Original EOS
                forceSCx=0.;
                forceSCy=0.;
                forceSCz=0.;
                for(pop=0;pop<170;pop++){
                    tempID2=((k+cz10[pop])*sizeX2*sizeY2)+((j+cy10[pop])*sizeX2)+(i+cx10[pop]);
                    rhoTemp2=rho[tempID2];
                    psiSC2=rhoSC*(1.-exp(0.-rhoTemp2/rhoSC)); //fluid-solid interaction strength determined by rhowall
                    forceSCx=forceSCx-(cs*cs*potGSC*weight10[pop]*psiSC1*psiSC2*(double)cx10[pop]);
                    forceSCy=forceSCy-(cs*cs*potGSC*weight10[pop]*psiSC1*psiSC2*(double)cy10[pop]);
                    forceSCz=forceSCz-(cs*cs*potGSC*weight10[pop]*psiSC1*psiSC2*(double)cz10[pop]);
                }
                //sum interfacial and exterior forces
                forceAllx[tempID]=forceSCx+accell[0]*rhoTemp1;
                forceAlly[tempID]=forceSCy+accell[1]*rhoTemp1;
                forceAllz[tempID]=forceSCz+accell[2]*rhoTemp1;
                //if(rank==14 && i==50 && k==1)printf("iT=%d,i=%d,j=%d,k=%d forces: %E %E %E psi1 = %E, psi2e18 = %E\n",iT,i,j,k,forceAll[0][i][j][k],forceAll[1][i][j][k],forceAll[2][i][j][k], psiSC1, psiSC2);
            }
        }


        ///Macroscopic Variables calculated on all fluid nodes with velocity shift
        for(iD=0; iD<=countF; iD++){
            tempID=listFnodes[iD];
            cxfIn=0;
            cyfIn=0;
            czfIn=0;
            rhoTempInv=1./rho[tempID];
            for(pop=0;pop<=18;pop++){                   //calculate velocity dot products and sum pops for density
                tempIDpop=sizepop*tempID+pop;
                cxfIn=cxfIn+cxD[pop]*fIn[tempIDpop];
                cyfIn=cyfIn+cyD[pop]*fIn[tempIDpop];
                czfIn=czfIn+czD[pop]*fIn[tempIDpop];
            }
            ux[tempID]=(cxfIn+0.5*forceAllx[tempID])*rhoTempInv; // forcing with Guo velocity shift
            uy[tempID]=(cyfIn+0.5*forceAlly[tempID])*rhoTempInv;
            uz[tempID]=(czfIn+0.5*forceAllz[tempID])*rhoTempInv;
        }

        ///fluid nodes undergo collision
        if(globalVisc_==1){
            ///global viscosity collision
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                jx=rho[tempID]*ux[tempID];
                jy=rho[tempID]*uy[tempID];
                jz=rho[tempID]*uz[tempID];
                rhoTemp=rho[tempID];
                mEq[0]=rhoTemp;
                mEq[1]=(-11.*rhoTemp)+((19./rhoFluid)*(jx*jx+jy*jy+jz*jz)); //rhoFluid is user input or calculated as mean of density if pressure gradient applied
                mEq[2]=(omEp*rhoTemp)+((omEpj/rhoFluid)*(jx*jx+jy*jy+jz*jz));
                mEq[3]=jx;
                mEq[4]=(-2./3.)*jx;
                mEq[5]=jy;
                mEq[6]=(-2./3.)*jy;
                mEq[7]=jz;
                mEq[8]=(-2./3.)*jz;
                mEq[9]=(1./rhoFluid)*((2.*jx*jx)-(jy*jy)-(jz*jz));
                mEq[10]=omxx*mEq[9];
                mEq[11]=(1./rhoFluid)*(jy*jy-jz*jz);
                mEq[12]=omxx*mEq[11];
                mEq[13]=(1./rhoFluid)*jx*jy;
                mEq[14]=(1./rhoFluid)*jy*jz;
                mEq[15]=(1./rhoFluid)*jz*jx;
                mEq[16]=0.;
                mEq[17]=0.;
                mEq[18]=0.;
                //force
                uFx=ux[tempID]*forceAllx[tempID];
                uFy=uy[tempID]*forceAlly[tempID];
                uFz=uz[tempID]*forceAllz[tempID];
                Fx=forceAllx[tempID];
                Fy=forceAlly[tempID];
                Fz=forceAllz[tempID];
                MF[0]=0.;
                MF[1]=38.*(uFx+uFy+uFz);
                MF[2]=-11.*(uFx+uFy+uFz);
                MF[3]=Fx;
                MF[4]=(-2./3.)*Fx;
                MF[5]=Fy;
                MF[6]=(-2./3.)*Fy;
                MF[7]=Fz;
                MF[8]=(-2./3.)*Fz;
                MF[9]=2.*(2.*uFx-uFy-uFz);
                MF[10]=-2.*uFx+uFy+uFz;
                MF[11]=2.*(uFy-uFz);
                MF[12]=-uFy+uFz;
                MF[13]=uy[tempID]*Fx+ux[tempID]*Fy;
                MF[14]=uz[tempID]*Fy+uy[tempID]*Fz;
                MF[15]=uz[tempID]*Fx+ux[tempID]*Fz;
                MF[16]=0.;
                MF[17]=0.;
                MF[18]=0.;
                //Collision
                for(popM=0;popM<=18;popM++){
                    m[popM]=0;
                    for(pop=0;pop<=18;pop++){ //transform distribution function (velocity space) to moment-space
                        tempIDpop=(sizepop*tempID)+pop;
                        m[popM] = m[popM] + M[popM][pop]*fIn[tempIDpop];
                    }
                    colM[popM] = S[popM]*(m[popM] - mEq[popM]); //collide with relaxation of source
                    colMF[popM] = (1.0-0.5*S[popM])*(MF[popM]);
                }
                for(pop=0;pop<=18;pop++){
                    f[pop]=0;
                    source[pop]=0.;
                    for(popM=0;popM<=18;popM++){ //transform moments back to velocity-space
                        f[pop] = f[pop] + invM[pop][popM]*colM[popM];
                        source[pop]=source[pop] + invM[pop][popM]*colMF[popM];
                    }
                    tempIDpop=(sizepop*tempID)+pop;
                    fOut[tempIDpop] = fIn[tempIDpop] - f[pop] + source[pop];
                }
            }
        }
        ///LEV model collision
        else{
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                jx=rho[tempID]*ux[tempID];
                jy=rho[tempID]*uy[tempID];
                jz=rho[tempID]*uz[tempID];
                rhoTemp=rho[tempID];
                S[9]=svLocSub[tempID]; ///This is LEV
                S[11]=S[9];
                S[13]=S[9];
                S[14]=S[9];
                S[15]=S[9];
                S[4]=sqLocSub[tempID];
                S[6]=S[4];
                S[8]=S[4];
                mEq[0]=rhoTemp;
                mEq[1]=(-11.*rhoTemp)+((19./rhoFluid)*(jx*jx+jy*jy+jz*jz)); //rhoFluid is user input or calculated as mean of density if pressure gradient applied
                mEq[2]=(omEp*rhoTemp)+((omEpj/rhoFluid)*(jx*jx+jy*jy+jz*jz));
                mEq[3]=jx;
                mEq[4]=(-2./3.)*jx;
                mEq[5]=jy;
                mEq[6]=(-2./3.)*jy;
                mEq[7]=jz;
                mEq[8]=(-2./3.)*jz;
                mEq[9]=(1./rhoFluid)*((2.*jx*jx)-(jy*jy)-(jz*jz));
                mEq[10]=omxx*mEq[9];
                mEq[11]=(1./rhoFluid)*(jy*jy-jz*jz);
                mEq[12]=omxx*mEq[11];
                mEq[13]=(1./rhoFluid)*jx*jy;
                mEq[14]=(1./rhoFluid)*jy*jz;
                mEq[15]=(1./rhoFluid)*jz*jx;
                mEq[16]=0.;
                mEq[17]=0.;
                mEq[18]=0.;
                //force
                uFx=ux[tempID]*forceAllx[tempID];
                uFy=uy[tempID]*forceAlly[tempID];
                uFz=uz[tempID]*forceAllz[tempID];
                Fx=forceAllx[tempID];
                Fy=forceAlly[tempID];
                Fz=forceAllz[tempID];
                MF[0]=0.;
                MF[1]=38.*(uFx+uFy+uFz);
                MF[2]=-11.*(uFx+uFy+uFz);
                MF[3]=Fx;
                MF[4]=(-2./3.)*Fx;
                MF[5]=Fy;
                MF[6]=(-2./3.)*Fy;
                MF[7]=Fz;
                MF[8]=(-2./3.)*Fz;
                MF[9]=2.*(2.*uFx-uFy-uFz);
                MF[10]=-2.*uFx+uFy+uFz;
                MF[11]=2.*(uFy-uFz);
                MF[12]=-uFy+uFz;
                MF[13]=uy[tempID]*Fx+ux[tempID]*Fy;
                MF[14]=uz[tempID]*Fy+uy[tempID]*Fz;
                MF[15]=uz[tempID]*Fx+ux[tempID]*Fz;
                MF[16]=0.;
                MF[17]=0.;
                MF[18]=0.;
                //Collision
                for(popM=0;popM<=18;popM++){
                    m[popM]=0;
                    for(pop=0;pop<=18;pop++){ //transform distribution function (velocity space) to moment-space
                        tempIDpop=(sizepop*tempID)+pop;
                        m[popM] = m[popM] + M[popM][pop]*fIn[tempIDpop];
                    }
                    colM[popM] = S[popM]*(m[popM] - mEq[popM]); //collide with relaxation of source
                    colMF[popM] = (1.0-0.5*S[popM])*(MF[popM]);
                }
                for(pop=0;pop<=18;pop++){
                    f[pop]=0;
                    source[pop]=0.;
                    for(popM=0;popM<=18;popM++){ //transform moments back to velocity-space
                        f[pop] = f[pop] + invM[pop][popM]*colM[popM];
                        source[pop]=source[pop] + invM[pop][popM]*colMF[popM];
                    }
                    tempIDpop=(sizepop*tempID)+pop;
                    fOut[tempIDpop] = fIn[tempIDpop] - f[pop] + source[pop];
                }
            }
        }


///++++++++++++++++++++++++++++CHECK FOR CONVERGENCE+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

        //Determine uDarcy for each processor and send to rank 0 to determine overall mean darcy velocity and check for convergence
        if (iT%tConv==0){
            rho_sum=0.;
            rhoAll=0.;
            uDarcy_sum=0;
            uAll=0;
            for(iD=0; iD<=countF; iD++){
                tempID=listFnodes[iD];
                if(flowDir==0){
                    uAll=uAll+ux[tempID];
                    rhoAll=rhoAll+rho[tempID];
                }
                else if (flowDir==1){
                    uAll=uAll+uy[tempID];
                    rhoAll=rhoAll+rho[tempID];
                }
                else{		//flowDir=2, flow in z-direction is the default
                    uAll=uAll+uz[tempID];
                    rhoAll=rhoAll+rho[tempID];
                }
            }
            if (rank>0){
                MPI_Send(&uAll, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);  //send all u and rho to rank 0
                MPI_Send(&rhoAll, 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
            }
            else{  //rank 0 collects all other processors uDarcy
                uDarcy_sum=uDarcy_sum+uAll;	//uAll
                rho_sum=rho_sum+rhoAll;
                //printf("uDarcy sum for rank %i:  %lf\n", rank, uAll);
                for (thread=1;thread<=size-1;thread++){
                    MPI_Recv(&uAll, 1, MPI_DOUBLE, thread, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&rhoAll, 1, MPI_DOUBLE, thread, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    uDarcy_sum=uDarcy_sum+uAll;
                    rho_sum=rho_sum+rhoAll;
                }
                //printf("uDarcy sum for all:  %lf\n", uDarcy_sum);
                uDarcy_mean=uDarcy_sum/((double)lx*(double)ly*(double)lz);
                uFlow_mean=uDarcy_sum/((double)countF_sum);
                rho_mean=rho_sum/((double)countF_sum);
                perm=uDarcy_mean*mu/gradP_perm;
                permFlow=uFlow_mean*mu/gradP_perm;
                walltime2=MPI_Wtime();
                printf("iteration number: %i\n", iT);
                printf("Simulation Run Time: %f\n", walltime2-walltime1);
                printf("Mean Darcy Velocity: %E\n", uDarcy_mean);
                printf("Mean Flow Velocity: %E\n", uFlow_mean);
                printf("Mean Density: %E\n", rho_mean);
                printf("Permeability (Darcy): %E\n", perm);
                printf("Permeability Pore Space Only: %E\n", permFlow);
                convU=sqrt((uDarcy_mean-uDarcy_old)*(uDarcy_mean-uDarcy_old))/sqrt(uDarcy_mean*uDarcy_mean);
                printf("Velocity convergence: %E\n", convU);
                if (convU < epsilon){
                    printf("simulation has reached velocity convergence\n");
                    //set iT to maxT to break while loop
                    iT=maxT;
                }
                //tell all other ranks to break if iT is set to maxT
                for(thread=1;thread<=size-1;thread++){
                    MPI_Send(&iT, 1, MPI_INT, thread, 4, MPI_COMM_WORLD);
                }
                uDarcy_old=uDarcy_mean;
            }
            if(rank>0){
                MPI_Recv(&iT, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
///++++++++++++++++++++++++++++ OUTPUT MACROSCALE OR MICROSCALE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

               //outputs
        if(iT%tOut==0 || iT>previT){
            walltime2=MPI_Wtime();
            if (rank==0){
                printf("\nAt iteration: %i \nSimulation Time [s]: %lf\n", previT, walltime2-walltime1);
            }
            mkdir(nameFile,0777);
            //output velocities
            if(outux_==1){
                char imageFilename1[150];
                char folderName1[150];
                sprintf(folderName1, "%s/%s_ux", nameFile, nameFile);
                mkdir(folderName1,0777);
                sprintf(imageFilename1, "%s/%s_ux/ux_%i_%i",nameFile,nameFile,rank,previT);
                FILE *imagef1 = fopen(imageFilename1,"w");
                for(k=isoNum;k<sizeZ2-isoNum;k++){
                    for(j=isoNum;j<sizeY2-isoNum;j++){
                        for(i=isoNum;i<sizeX2-isoNum;i++){
                            tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                            outTemp = ux[tempID];
                            fprintf(imagef1, "%.18f\n", outTemp);
                        }
                    }
                }
                fclose(imagef1);
            }
            if(outuy_==1){
                char imageFilename2[150];
                char folderName2[150];
                sprintf(folderName2, "%s/%s_uy", nameFile, nameFile);
                mkdir(folderName2,0777);
                sprintf(imageFilename2, "%s/%s_uy/uy_%i_%i",nameFile,nameFile,rank,previT);
                FILE *imagef2 = fopen(imageFilename2,"w");
                for(k=isoNum;k<sizeZ2-isoNum;k++){
                    for(j=isoNum;j<sizeY2-isoNum;j++){
                        for(i=isoNum;i<sizeX2-isoNum;i++){
                            tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                            outTemp = uy[tempID];
                            fprintf(imagef2, "%.18f\n", outTemp);
                        }
                    }
                }
                fclose(imagef2);
            }
            if(outuz_==1){
                char imageFilename3[150];
                char folderName3[150];
                sprintf(folderName3, "%s/%s_uz", nameFile, nameFile);
                mkdir(folderName3,0777);
                sprintf(imageFilename3, "%s/%s_uz/uz_%i_%i",nameFile,nameFile,rank,previT);
                FILE *imagef3 = fopen(imageFilename3,"w");
                for(k=isoNum;k<sizeZ2-isoNum;k++){
                    for(j=isoNum;j<sizeY2-isoNum;j++){
                        for(i=isoNum;i<sizeX2-isoNum;i++){
                            tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                            outTemp = uz[tempID];
                            fprintf(imagef3, "%.18f\n", outTemp);
                        }
                    }
                }
                fclose(imagef3);
            }

            //output densities
            if(outrho_==1){
                char imageFilename4[150];
                char folderName4[150];
                sprintf(folderName4, "%s/%s_rho", nameFile,nameFile);
                mkdir(folderName4,0777);
                sprintf(imageFilename4, "%s/%s_rho/rho_%i_%i",nameFile,nameFile,rank,previT);
                FILE *imagef4 = fopen(imageFilename4,"w");
                for(k=isoNum;k<sizeZ2-isoNum;k++){
                    for(j=isoNum;j<sizeY2-isoNum;j++){
                        for(i=isoNum;i<sizeX2-isoNum;i++){
                            tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                            outTemp = rho[tempID];
                            fprintf(imagef4, "%.18f\n", outTemp);
                        }
                    }
                }
                fclose(imagef4);
            }

            //output particle distributions
            if(outfIn_==1){
                for(pop=0;pop<=18;pop++){
                    char imageFilename5[150];
                    char folderName5[150];
                    sprintf(folderName5, "%s/%s_fIn", nameFile, nameFile);
                    mkdir(folderName5,0777);
                    sprintf(imageFilename5, "%s/%s_fIn/fIn_%i_%i_%i",nameFile, nameFile,rank,pop,previT);
                    FILE *imagef5 = fopen(imageFilename5,"w");
                    for(k=isoNum;k<sizeZ2-isoNum;k++){
                        for(j=isoNum;j<sizeY2-isoNum;j++){
                            for(i=isoNum;i<sizeX2-isoNum;i++){
                                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                                outTemp = fIn[tempIDpop];
                                fprintf(imagef5, "%.18f\n", outTemp);
                            }
                        }
                    }
                    fclose(imagef5);
                }
            }
            if(outfOut_==1){
                for(pop=0;pop<=18;pop++){
                    char imageFilename6[150];
                    char folderName6[150];
                    sprintf(folderName6, "%s/%s_fOut", nameFile,nameFile);
                    mkdir(folderName6,0777);
                    sprintf(imageFilename6, "%s/%s_fOut/fOut_%i_%i_%i",nameFile,nameFile,rank,pop,previT);
                    FILE *imagef6 = fopen(imageFilename6,"w");
                    for(k=isoNum;k<sizeZ2-isoNum;k++){
                        for(j=isoNum;j<sizeY2-isoNum;j++){
                            for(i=isoNum;i<sizeX2-isoNum;i++){
                                tempIDpop=(k*sizepop*sizeX2*sizeY2)+(j*sizepop*sizeX2)+(i*sizepop)+pop;
                                outTemp = fOut[tempIDpop];
                                fprintf(imagef6, "%.18f\n", outTemp);
                            }
                        }
                    }
                    fclose(imagef6);
                }
            }

            //output geometry
            if(outgeom_==1){
                char imageFilename7[150];
                char folderName7[150];
                sprintf(folderName7, "%s/%s_GeometryOut", nameFile, nameFile);
                mkdir(folderName7,0777);
                sprintf(imageFilename7, "%s/%s_GeometryOut/GeometryOut_%i",nameFile,nameFile,rank);
                FILE *imagef7 = fopen(imageFilename7,"w");
                for(k=isoNum;k<sizeZ2-isoNum;k++){
                    for(j=isoNum;j<sizeY2-isoNum;j++){
                        for(i=isoNum;i<sizeX2-isoNum;i++){
                            tempID=(k*sizeX2*sizeY2)+(j*sizeX2)+i;
                            outTempInt = bbRegionSub[tempID];
                            fprintf(imagef7, "%d\n", outTempInt);
                        }
                    }
                }
            fclose(imagef7);
            }
        }

    }

    MPI_Finalize();
    return 0;
}

