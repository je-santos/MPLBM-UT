#include "palabos3D.h"
#include "palabos3D.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace plb;
using namespace std;

typedef double T; // Use double-precision arithmetics
#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor // Use a grid which additionally to the f's stores two variables for the external force term.

void writeGifs(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,  //creates the pictures
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, plint runs, plint iT)
{
    const plint imSize = 600;
    const plint zcomponent = 0;
    const plint nx = lattice_fluid2.getNx();
    const plint ny = lattice_fluid2.getNy();
    const plint nz = lattice_fluid2.getNz();
    Box3D slice(0, nx, 0, ny, nz / 2, nz / 2);

    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("rho_f1_", 100000000 * runs + iT, 8),
        *computeDensity(lattice_fluid1, slice), imSize, imSize);
    //imageWriter.writeScaledGif(createFileName("rho_f2_", iT, 8),*computeDensity(lattice_fluid2, slice), imSize,imSize);
    //imageWriter.writeScaledGif(createFileName("velz_f1", iT, 8),*computeVelocityComponent(lattice_fluid1, slice, zcomponent), imSize,imSize);
}


void writeVTK_vel(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, plint runs)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk_vel_rho1_", 1 * runs, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice_fluid1), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice_fluid1), "velocity", 1.);
}



void writeGifs2(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, plint runs, plint iT)
{
    const plint imSize = 600;
    const plint zcomponent = 0;
    const plint nx = lattice_fluid2.getNx();
    const plint ny = lattice_fluid2.getNy();
    const plint nz = lattice_fluid2.getNz();
    Box3D slice(0, nx, ny / 2, ny / 2, 0, nz);

    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("y_rho_f1_", 1 * runs, 8),
        *computeDensity(lattice_fluid1, slice), imSize, imSize);
    //imageWriter.writeScaledGif(createFileName("rho_f2_", iT, 8),*computeDensity(lattice_fluid2, slice), imSize,imSize);
    //imageWriter.writeScaledGif(createFileName("velz_f1", iT, 8),*computeVelocityComponent(lattice_fluid1, slice, zcomponent), imSize,imSize);
}

void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, plint runs, plint iter, plint nx, plint ny, plint nz) // std::string Dstring, std::string Vstring
{
    const plint zcomponent = 0;

    VtkImageOutput3D<double> vtkOut(createFileName("rho1_", 1* runs, 8), 1.);
    vtkOut.writeData<double>((*computeDensity(lattice)), "Density", 1.);

    //VtkImageOutput3D<double> vtkvel_Out(createFileName("vtk_vel", 100000000*runs+iter, 8), 1.);
    //vtkvel_Out.writeData<double>((*computeVelocityComponent(lattice,Box3D(0, nx-1,0,ny-1,0,nz-1),zcomponent)), "z velocity", 1.);
}

void writeVTK2(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, plint runs, plint iter, plint nx, plint ny, plint nz) // std::string Dstring, std::string Vstring
{
    const plint xcomponent = 0;

    VtkImageOutput3D<double> vtkOut(createFileName("rho2_", 100000000 * runs + iter, 8), 1.);
    vtkOut.writeData<double>((*computeDensity(lattice)), "Density", 1.);

    //VtkImageOutput3D<double> vtkvel_Out(createFileName("vtk_vel", iter, 8), 1.);
    //vtkvel_Out.writeData<double>((*computeVelocityComponent(lattice,Box3D(0, nx-1,0,ny-1,0,nz-1),xcomponent)), "x velocity", 1.);
}

T computeVelocity1(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, T nu1)
{
	plint xComponent = 0;
	const plint nx = lattice_fluid1.getNx();
    const plint ny = lattice_fluid1.getNy();
    const plint nz = lattice_fluid1.getNz();
	
    Box3D domain(3, nx-3, 0, ny-1, 0, nz-1);
	T meanU1 = computeAverage(*computeVelocityComponent(lattice_fluid1, domain, xComponent));
	pcout << "Average velocity for fluid1 in x direction    = " << meanU1            << std::endl;
	return meanU1;
}

T computeVelocity2(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid2, T nu2)
{
	plint xComponent = 0;
	const plint nx = lattice_fluid2.getNx();
    const plint ny = lattice_fluid2.getNy();
    const plint nz = lattice_fluid2.getNz();
	
    Box3D domain(3, nx-3, 0, ny-1, 0, nz-1);
	T meanU2 = computeAverage(*computeVelocityComponent(lattice_fluid2, domain, xComponent));
	pcout << "Average velocity for fluid2 in x direction    = " << meanU2            << std::endl;
	return meanU2;
}


void PorousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, MultiScalarField3D<int>& geometry,
    T rhoNoFluid, T rho_fluid1, T rho_fluid2, T Gads_f1_s1, T Gads_f1_s2, T force_fluid1, T force_fluid2, T nx1f1, T nx2f1, T ny1f1, T ny2f1, T nz1f1, T nz2f1, T nx1f2, T nx2f2, T ny1f2, T ny2f2, T nz1f2, T nz2f2)
{
    plint nx = lattice_fluid2.getNx();
    plint ny = lattice_fluid2.getNy();
    plint nz = lattice_fluid2.getNz();

    pcout << "Definition of the geometry." << endl;
    for (plint itX = 0; itX < nx; ++itX) {
        for (plint itY = 0; itY < ny; ++itY) {
            for (plint itZ = 0; itZ < nz; ++itZ) {
                if (geometry.get(itX, itY, itZ) == 1) { //wall 1
                    defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(Gads_f1_s1));
                    defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-Gads_f1_s1));
                }
                if (geometry.get(itX, itY, itZ) == 2) { //solid
                    defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
                    defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
                }
                if (geometry.get(itX, itY, itZ) == 3) { //wall 2
                    defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(Gads_f1_s2));
                    defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-Gads_f1_s2));
                }
                if (geometry.get(itX, itY, itZ) == 4) {
                    defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-0.3));
                    defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(0.3));
                }
                if (geometry.get(itX, itY, itZ) == 5) {
                    defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(0));
                    defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(0));
                }
            }
        }
    }

    // Output geometry dynamics
    VtkImageOutput3D<int> vtkOut(createFileName("vtkgeometry", 1, 1), 1.);
    vtkOut.writeData<int>(geometry, "Dynamics", 1.);

    Array<T, 3> zeroVelocity(0., 0., 0.);

    // Initialize  uniform density for target saturation
    pcout << "Initializing Fluids" << endl;

    initializeAtEquilibrium(lattice_fluid2, Box3D(nx1f2, nx2f2, ny1f2, ny2f2, nz1f2, nz2f2), rho_fluid2, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid1, Box3D(nx1f2, nx2f2, ny1f2, ny2f2, nz1f2, nz2f2), rhoNoFluid, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid1, Box3D(nx1f1, nx2f1, ny1f1, ny2f1, nz1f1, nz2f1), rho_fluid1, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid2, Box3D(nx1f1, nx2f1, ny1f1, ny2f1, nz1f1, nz2f1), rhoNoFluid, zeroVelocity);

    //initializeAtEquilibrium(lattice_fluid1, Box3D(0, 4, 0, ny, 0, nz), rho_fluid1, zeroVelocity);
    //initializeAtEquilibrium(lattice_fluid2, Box3D(0, 4, 0, ny, 0, nz), rhoNoFluid, zeroVelocity);

    //initializeAtEquilibrium(lattice_fluid2, Box3D(5, nx, 0, ny, 0, nz), rho_fluid2, zeroVelocity);
    //initializeAtEquilibrium(lattice_fluid1, Box3D(5, nx, 0, ny, 0, nz), rhoNoFluid, zeroVelocity);

    setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_fluid1, 0.));
    setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_fluid2, 0.));

    lattice_fluid1.initialize();
    lattice_fluid2.initialize();
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    std::string fNameOut = argv[1];
    global::directories().setOutputDir(fNameOut);
    const T omega1 = atof(argv[2]);
    const T omega2 = atof(argv[3]);
    const plint nx = atoi(argv[4]);
    const plint ny = atoi(argv[5]);
    const plint nz = atoi(argv[6]);
    const T G = atof(argv[7]);
    const T force_fluid1 = atof(argv[8]);
    const T force_fluid2 = atof(argv[9]);
    const T Gads_f1_s1 = atof(argv[10]);
    const T Gads_f1_s2 = atof(argv[11]);
    std::string fNameIn = argv[12];
    const plint maxIter = atoi(argv[13]); //max_no_it
    const plint saveIter = atoi(argv[14]);
    const plint vtkIter = atoi(argv[15]);
    const plint statIter = atoi(argv[16]);
    const T convergence = atof(argv[17]);
    const T tempHold = 1000000;
    const plint startNum = atoi(argv[18]);
    //const plint runnum = atoi(argv[19]);// no of density simulations to run
    const T nx1f1 = atoi(argv[20]); //fluid configuration
    const T nx2f1 = atoi(argv[21]);
    const T ny1f1 = atoi(argv[22]);
    const T ny2f1 = atoi(argv[23]);
    const T nz1f1 = atoi(argv[24]);
    const T nz2f1 = atoi(argv[25]);
    const T nx1f2 = atoi(argv[26]);
    const T nx2f2 = atoi(argv[27]);
    const T ny1f2 = atoi(argv[28]);
    const T ny2f2 = atoi(argv[29]);
    const T nz1f2 = atoi(argv[30]);
    const T nz2f2 = atoi(argv[31]);
    const T rhoNoFluid = atof(argv[32]);
    const T cap_factor = atof(argv[33]);
	const T rho_fluid1_min = atof(argv[34]);
	const T rho_fluid1_max = atof(argv[35]);
	const T rho_fluid1_step = atof(argv[36]);
	const T rho_fluid2_val = atof(argv[37]);
	const plint runnum = ((rho_fluid1_max - rho_fluid1_min) / rho_fluid1_step)+1;	
    T rho_fluid1[runnum] ;
    T rho_fluid2[runnum];
	T M[runnum];
	T Ca[runnum];
	T Mom1[runnum];
	T Mom2[runnum];
	T Kr1[runnum];
	T Kr2[runnum];
	T Mom1_high;
	T Mom2_high;
	T deltaP[runnum];
	T kr1[runnum];
	T kr2[runnum];
	T k1_high;
	T k2_high;
	/*
    pcout<<"rho_fluid1"<<endl;
	for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid1[readnum] = atof(argv[(34 + ((readnum - 1)))]);
		pcout << rho_fluid1[readnum];
    }
	pcout<<endl;
	pcout<<"rho_fluid2"<<endl;
    for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid2[readnum] = atof(argv[(34 + runnum + ((readnum - 1)))]);
		pcout << rho_fluid2[readnum];
	}
	pcout<<endl;
	*/
	
	
	
	for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid1[readnum] = rho_fluid1_min + (readnum-1)* rho_fluid1_step;
    }
	for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid2[readnum] = rho_fluid2_val;
	}
	
    const T nu1 = ((T)1 / omega1 - 0.5) / DESCRIPTOR<T>::invCs2;
    const T nu2 = ((T)1 / omega2 - 0.5) / DESCRIPTOR<T>::invCs2;

    // Use regularized BGK dynamics to improve numerical stability (but note that BGK dynamics works well too).
    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega2));
    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1));

    lattice_fluid2.periodicity().toggle(0, false);
    lattice_fluid1.periodicity().toggle(0, false);
    lattice_fluid2.periodicity().toggle(1, false);
    lattice_fluid1.periodicity().toggle(1, false);
    lattice_fluid2.periodicity().toggle(2, false);
    lattice_fluid1.periodicity().toggle(2, false);

    // Store a pointer to all lattices (two in the present application) in a vector to
    //   create the Shan/Chen coupling term. The fluid 2 being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to collideAndStream() or stream() for the fluid 2.
    vector<MultiBlockLattice3D<T, DESCRIPTOR>*> blockLattices;
    blockLattices.push_back(&lattice_fluid2);
    blockLattices.push_back(&lattice_fluid1);

    // The argument "constOmegaValues" to the Shan/Chen processor is optional,
    //   and is used for efficiency reasons only. It tells the data processor
    //   that the relaxation times are constant, and that their inverse must be
    //   computed only once.
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega2);
    constOmegaValues.push_back(omega1);
    plint processorLevel = 1;
    integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D<T, DESCRIPTOR>(G, constOmegaValues), Box3D(0, nx, 0, ny, 0, nz), blockLattices, processorLevel); // efficiency reasons

    pcout << "Convergence = " << convergence << endl;
    pcout << "nx = " << nx << endl;
    pcout << "ny = " << ny << endl;
    pcout << "nz = " << nz << endl;
    pcout << "Gc = " << G << endl;
    pcout << "force fluid 1 = " << force_fluid1 << endl;
    pcout << "force fluid 2 = " << force_fluid2 << endl;
    pcout << "G_ads_1 fluid 1 = " << Gads_f1_s1 << endl;
    pcout << "G_ads_1 fluid 2 = " << -Gads_f1_s1 << endl;
    pcout << "G_ads_2 fluid 1 = " << Gads_f1_s2 << endl;
    pcout << "G_ads_2 fluid 2 = " << -Gads_f1_s2 << endl;
    pcout << "nz2f2  = " << nz2f2 << endl;
    pcout << "Rho_no_fluid = " << rhoNoFluid << endl;

	    for (plint readnum = 1; readnum <= runnum; ++readnum) {
        pcout << "Rho_no_1 = " << rho_fluid1[readnum] << endl;
        pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
    }

    pcout << "Cap Multiplier = " << cap_factor << endl;
    pcout << "Reading the geometry file." << endl;

    MultiScalarField3D<int> geometry(nx, ny, nz);
    plb_ifstream geometryFile(fNameIn.c_str());

    if (!geometryFile.is_open()) {
        pcout << "Error: could not open geometry file " << fNameIn << endl;
        return -1;
    }
    geometryFile >> geometry;

    util::ValueTracer<T> converge1(1.0, lattice_fluid2.getNx(), convergence); // convergence parameters velocity/size/threshold
    util::ValueTracer<T> converge2(1.0, lattice_fluid2.getNx(), convergence);

    for (plint runs = 1; runs <= runnum; ++runs) { // Loop simulations with varying saturation

        if (runs > 1) {
            pcout << "Using previous Geometry  " << endl;
        }
        else {
            PorousMediaSetup(lattice_fluid1, lattice_fluid2, geometry, rhoNoFluid, rho_fluid1[1], rho_fluid2[1], Gads_f1_s1, Gads_f1_s2, force_fluid1, force_fluid2, nx1f1, nx2f1, ny1f1, ny2f1, nz1f1, nz2f1, nx1f2, nx2f2, ny1f2, ny2f2, nz1f2, nz2f2);
        }

        T meanJ_2 = 10000.0;
        T meanJ_old2;
        T meanJ_1 = 10000.0;
        T meanJ_old1;
        T convergemeanJ1;
        T convergemeanJ2;
        T convergedmeanJ1;
        T convergedmeanJ2;
		// Calculate dynamic viscosity ratio here and pcout 
		const T mu1 = rho_fluid1[runs]*nu1;
		const T mu2 = rho_fluid2[runs]*nu2;
		T M = mu1 / mu2;
		
        pcout << endl
              << "Starting simulation with rho 1:  " << rho_fluid1[runs] << endl;
        pcout << endl
              << "Starting simulation with rho 2:  " << rho_fluid2[runs] << endl;

        plint checkconv = 0;
        plint iT = 0;
        while (checkconv == 0) { // Main loop over time iterations.
            iT = iT + 1;

            // Time iteration for the fluid 1.
            lattice_fluid1.collideAndStream();
            // Time iteration for the fluid 2 must come after the fluid 1,
            //   because the coupling is executed here. You should understand this as follows.
            //   The effect of the coupling is to compute the interaction force between
            //   species, and to precompute density and momentum for each species. This must
            //   be executed *before* collide-and-streaming the fluids, because the collision
            //   step needs to access all these values. In the present case, it is done after
            //   both collide-and-stream step, which means, before the collide-and-stream of
            //   the next iteration (it's the same if you are before or after; the important
            //   point is not to be between the two collide-and-streams of the light and heavy
            //   fluid. As for the initial condition, the coupling is initially performed once
            //   during the function call to lattice_fluid2.initialize().
            lattice_fluid2.collideAndStream();

            if (iT % saveIter == 0) {
                writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
            }
            if (iT % vtkIter == 0) {
                writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
				writeVTK_vel(lattice_fluid1, runs + startNum);
                writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
				
            }

            if (iT == 1) {
                //writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                //writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                //writeVTK(lattice_fluid1, runs + startNum,iT, nx, ny, nz);
                //writeVTK2(lattice_fluid2, runs + startNum,iT, nx, ny, nz);
            }

converge1.takeValue(getStoredAverageEnergy(lattice_fluid1), true); //check for convergence
converge2.takeValue(getStoredAverageEnergy(lattice_fluid2), true); //check for convergence

            if (iT % statIter == 0) {
                pcout << "Iteration:  " << iT << endl;
                //pcout << "Average density fluid one = " << getStoredAverageDensity<T>(lattice_fluid1) << endl;
                //pcout << "Average density fluid two = " << getStoredAverageDensity<T>(lattice_fluid2) << endl << endl;
                //pcout << "meanJ:fluid1" << endl;
                //meanJ_old1 = meanJ_1;
                //meanJ_1 = getStoredAverageEnergy<T>(lattice_fluid1);
                //pcout << meanJ_1 << endl;
                //convergemeanJ1 = std::sqrt((meanJ_1 - meanJ_old1) * (meanJ_1 - meanJ_old1));
                //convergedmeanJ1 = std::sqrt((convergence * meanJ_1) * (convergence * meanJ_1));
                //pcout << "Difference in Convergence criteria (sum J fluid1):  " << convergemeanJ1 << endl;
                //pcout << "Convergence criteria (sum J fluid1):  " << convergedmeanJ1 << endl;
                //pcout << "meanJ:fluid2" << endl;
                //meanJ_old2 = meanJ_2;
                //meanJ_2 = getStoredAverageEnergy<T>(lattice_fluid2);
                //pcout << meanJ_2 << endl;
                //convergemeanJ2 = std::sqrt((meanJ_2 - meanJ_old2) * (meanJ_2 - meanJ_old2));
                //convergedmeanJ2 = std::sqrt((convergence * meanJ_2) * (convergence * meanJ_2));
                //pcout << "Difference in Convergence criteria (sum J fluid2):  " << convergemeanJ2 << endl;
                //pcout << "Convergence criteria (sum J fluid2):  " << convergedmeanJ2 << endl << endl;

                if ((converge1.hasConverged()) && (converge2.hasConverged())) {
					pcout << "Converge 1 and Converge 2 converged for iteration:" << iT << endl; 
                    checkconv = 1;
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
					writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
					
				// Calculate velocity here for both fluids in x-direction and pcout
				computeVelocity1(lattice_fluid1, nu1);
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				
				// Calculate rel perm from Darcy law
				deltaP[runs]=(rho_fluid1[runs]-rho_fluid2[runs])/3;
				kr1[runs]= nu1*meanU1/(deltaP[runs]/(T)(nx-6));
				kr2[runs]= nu2*meanU2/(deltaP[runs]/(T)(nx-6));
				
				// calculate capillary number
				const T cosTheta = abs((Gads_f1_s2-Gads_f1_s1)/(G*(rho_fluid1[runs]-rhoNoFluid)*0.5));
				Ca[runs] = ((mu1*meanU1)/(0.15*cosTheta));
				
				// calculate rel perm from momentum for for both fluids (measure 1st and last step to divide and calculate rel perm)
				Mom1[runs]=rho_fluid1[runs]*meanU1;
				Mom2[runs]=rho_fluid2[runs]*meanU2;
				if (runs == 1) {
					Mom2_high = Mom2[runs];
					k2_high = kr2[runs];
				}
				if (runs == runnum) {
					Mom1_high = Mom1[runs];
					k1_high = kr1[runs];
				}
                }


                if ((convergemeanJ1 < convergedmeanJ1) && (convergemeanJ2 < convergedmeanJ2)) {
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    pcout << "Simulation has converged" << endl;
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
					
				// Calculate velocity here for both fluids in x-direction and pcout
				computeVelocity1(lattice_fluid1, nu1);
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				
				// Calculate rel perm from Darcy law
				deltaP[runs]=(rho_fluid1[runs]-rho_fluid2[runs])/3;
				kr1[runs]= nu1*meanU1/(deltaP[runs]/(T)(nx-6));
				kr2[runs]= nu2*meanU2/(deltaP[runs]/(T)(nx-6));
				
				// calculate capillary number
				const T cosTheta = abs((Gads_f1_s2-Gads_f1_s1)/(G*(rho_fluid1[runs]-rhoNoFluid)*0.5));
				Ca[runs] = ((mu1*meanU1)/(0.15*cosTheta));
				
				// calculate rel perm from momentum for for both fluids (measure 1st and last step to divide and calculate rel perm)
				Mom1[runs]=rho_fluid1[runs]*meanU1;
				Mom2[runs]=rho_fluid2[runs]*meanU2;
				if (runs == 1) {
					Mom2_high = Mom2[runs];
					k2_high = kr2[runs];
				}
				if (runs == runnum) {
					Mom1_high = Mom1[runs];
					k1_high = kr1[runs];
				}
                    //checkconv = 1;
                }
            }

            if (maxIter == iT) {
                pcout << "Simulation has reached maximum iteration" << endl;
                checkconv = 1;
                writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
                writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
                writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
				// Calculate velocity here for both fluids in x-direction and pcout
				computeVelocity1(lattice_fluid1, nu1);
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				
				// Calculate rel perm from Darcy law
				deltaP[runs]=(rho_fluid1[runs]-rho_fluid2[runs])/3;
				kr1[runs]= nu1*meanU1/(deltaP[runs]/(T)(nx-6));
				kr2[runs]= nu2*meanU2/(deltaP[runs]/(T)(nx-6));
				
				// calculate capillary number
				const T cosTheta = abs((Gads_f1_s2-Gads_f1_s1)/(G*(rho_fluid1[runs]-rhoNoFluid)*0.5));
				Ca[runs] = ((mu1*meanU1)/(0.15*cosTheta));
				
				// calculate rel perm from momentum for for both fluids (measure 1st and last step to divide and calculate rel perm)
				Mom1[runs]=rho_fluid1[runs]*meanU1;
				Mom2[runs]=rho_fluid2[runs]*meanU2;
				if (runs == 1) {
					Mom2_high = Mom2[runs];
					k2_high = kr2[runs];
				}
				if (runs == runnum) {
					Mom1_high = Mom1[runs];
					k1_high = kr1[runs];
				}
            }
            initializeAtEquilibrium(lattice_fluid1, Box3D(1, 2, 0, ny-1, 0, nz-1), rho_fluid1[runs] * cap_factor, Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid2, Box3D(1, 2, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid1, Box3D(nx - 3, nx-2, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid2, Box3D(nx - 3, nx-2, 0, ny-1, 0, nz-1), rho_fluid2[runs], Array<T, 3>(0., 0., 0.));

            lattice_fluid1.initialize();
            lattice_fluid2.initialize();
        }
    }
	// calculating relative permeability & outputting variables
	for (plint runs = 1; runs <= runnum; ++runs) {
	Kr1[runs] =	Mom1[runs] / Mom1_high;
	Kr2[runs] =	Mom2[runs] / Mom2_high;
	kr1[runs] = kr1[runs] / k1_high;
	kr2[runs] = kr2[runs] / k2_high;
	pcout << "Run    = " << runs            << std::endl;
	pcout << "Pressure difference =  " << deltaP[runs] << std::endl;
	pcout << "Viscosity ratio =  " << M[runs] << std::endl;
	pcout << "Capillary number =  " << Ca[runs] << std::endl;
	pcout << "Kr1 from momentum   = " << Kr1[runs]            << std::endl;	
	pcout << "Kr2 from momentum   = " << Kr2[runs]            << std::endl;
    pcout << "Kr1 from Darcy law   = " << kr1[runs]            << std::endl;	
	pcout << "Kr2 from Darcy law   = " << kr2[runs]            << std::endl;
	}
}
