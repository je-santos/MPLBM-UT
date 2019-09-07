// Modified to calculate velocity /density (momentum) from pressure difference (rho 2 decreasing), with single-phase flow included in runs 1 and 2 and correcting for spurious velocities in even runs -- working

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <time.h>

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

}

void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1, plint runs, plint iter, plint nx, plint ny, plint nz) // std::string Dstring, std::string Vstring
{
    const plint zcomponent = 0;

    VtkImageOutput3D<double> vtkOut(createFileName("rho1_", 1* runs, 8), 1.);
    vtkOut.writeData<double>((*computeDensity(lattice_fluid1)), "Density", 1.);


}

void writeVTK2(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, plint runs, plint iter, plint nx, plint ny, plint nz) // std::string Dstring, std::string Vstring
{
    const plint xcomponent = 0;

    VtkImageOutput3D<double> vtkOut(createFileName("rho2_", 100000000 * runs + iter, 8), 1.);
    vtkOut.writeData<double>((*computeDensity(lattice_fluid2)), "Density", 1.);


}

T computeVelocity1(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, T nu1)
{
	plint xComponent = 0;
	const plint nx = lattice_fluid1.getNx();
    const plint ny = lattice_fluid1.getNy();
    const plint nz = lattice_fluid1.getNz();
	
    Box3D domain(3, nx-4, 0, ny-1, 0, nz-1);
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
	
    Box3D domain(3, nx-4, 0, ny-1, 0, nz-1);
	T meanU2 = computeAverage(*computeVelocityComponent(lattice_fluid2, domain, xComponent));
	pcout << "Average velocity for fluid2 in x direction    = " << meanU2            << std::endl;
	return meanU2;
}

void computeRatios(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid2, T rho_F1, T rho_F2, T nu1, T nu2, T Gads_f1_s2, T Gads_f1_s1, T G, T rhoNoFluid, T meanU1, T& meanRho1, T& meanRho2, T& M, T& Ca)
{

	const plint nx = lattice_fluid1.getNx();
    const plint ny = lattice_fluid1.getNy();
    const plint nz = lattice_fluid1.getNz();
	
	Box3D domain(3, nx-4, 0, ny-1, 0, nz-1);
	
	
	meanRho1 = computeAverageDensity(lattice_fluid1, domain);
    meanRho2 = computeAverageDensity(lattice_fluid2, domain);
	

	// meanRho1 = rho_F1; // density of input and output
	// meanRho2 = rho_F2;	
	
	T mu1 = meanRho1*nu1;
	T mu2 = meanRho2*nu2;
	
	M = (mu1 / mu2);
	const T sigma = 0.15;			
	const T cosTheta = abs((Gads_f1_s2-Gads_f1_s1)/(G*(meanRho1-rhoNoFluid)*0.5));
	Ca = (abs(mu1*meanU1)/(sigma*cosTheta));
	
	meanRho1 = getStoredAverageDensity<T>(lattice_fluid1);
	meanRho2 = getStoredAverageDensity<T>(lattice_fluid2);

	pcout << "MeanRho1    = " << meanRho1            << std::endl;
	pcout << "MeanRho2    = " << meanRho2            << std::endl;
	pcout << "Domain 1     = " << nx <<"x" << ny <<"x" << nz            << std::endl;
	pcout << "Domain 2     = " << nx <<"x" << ny <<"x" << nz            << std::endl;	
	pcout << "M     = " << M            << std::endl;
	pcout << "Ca    = " << Ca            << std::endl;
}
void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>& geometry)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);
    std::auto_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
    plb_ifstream geometryFile(fNameIn.c_str());
    for (plint iX=0; iX<nx-1; ++iX) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << fNameIn << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryFile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
    }

    {
        VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
        vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
    }

    {
        std::auto_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int,T>(geometry, geometry.getBoundingBox());
        std::vector<T> isoLevels;
        isoLevels.push_back(0.5);
        typedef TriangleSet<T>::Triangle Triangle;
        std::vector<Triangle> triangles;
        Box3D domain = floatTags->getBoundingBox().enlarge(-1);
        domain.x0++;
        domain.x1--;
        isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
        TriangleSet<T> set(triangles);
        std::string outDir = fNameOut + "/";
        set.writeBinarySTL(outDir + "porousMedium.stl");
    }
}


 void PorousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
  MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, MultiScalarField3D<int>& geometry,
T rhoNoFluid, T rho1, T rho2, T Gads_f1_s1, T Gads_f1_s2, T force_fluid1, T force_fluid2, T nx1f1, T nx2f1, T ny1f1, T ny2f1, T nz1f1, T nz2f1, T nx1f2, T nx2f2, T ny1f2, T ny2f2, T nz1f2, T nz2f2, T runs, T new_file)
{
    plint nx = lattice_fluid2.getNx();
    plint ny = lattice_fluid2.getNy();
    plint nz = lattice_fluid2.getNz();
	const T rhoZero = 0.00001;
	const T deltaP = 0.05;

    pcout << "Definition of the geometry." << endl;
	if (new_file == 0) {	
	
			loadBinaryBlock(lattice_fluid1, "lattice_fluid1.dat");
			loadBinaryBlock(lattice_fluid2, "lattice_fluid2.dat");
	}
	
	if (new_file == 1) {
	
                defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>(Gads_f1_s1), 1);
                defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s1), 1);
				
                defineDynamics(lattice_fluid1, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);
                defineDynamics(lattice_fluid2, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

                defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>(Gads_f1_s2), 3);
                defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s2), 3);

                defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>(0), 4);
                defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(0), 4);

                defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>(0), 5);
                defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(0), 5);					
	
	
	
	}


    // Output geometry dynamics
	if (new_file == 1) {
    VtkImageOutput3D<int> vtkOut(createFileName("vtkgeometry", 1, 1), 1.);
    vtkOut.writeData<int>(geometry, "Dynamics", 1.);
	pcout << "Creating geometry vtk file" << endl;
	}

    Array<T, 3> zeroVelocity(0., 0., 0.);
			
	if (new_file == 1) {

    // Initialize  uniform density for target saturation
    pcout << "Initializing Fluids" << endl;

    initializeAtEquilibrium(lattice_fluid2, Box3D(nx1f2, nx2f2, ny1f2, ny2f2-1, nz1f2, nz2f2-1), rho2, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid1, Box3D(nx1f2, nx2f2, ny1f2, ny2f2-1, nz1f2, nz2f2-1), rhoNoFluid, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid1, Box3D(nx1f1-1, nx2f1, ny1f1, ny2f1-1, nz1f1, nz2f1-1), rho1, zeroVelocity);
    initializeAtEquilibrium(lattice_fluid2, Box3D(nx1f1-1, nx2f1, ny1f1, ny2f1-1, nz1f1, nz2f1-1), rhoNoFluid, zeroVelocity);
	
	 }


    setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_fluid1, 0.));
    setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_fluid2, 0.));

    lattice_fluid1.initialize();
    lattice_fluid2.initialize();
}

int main(int argc, char* argv[])
{
	clock_t t;
	t = clock();
    plbInit(&argc, &argv);
    // std::string fNameOut = argv[1];	
	  std::string fNameOut;
	  T new_file ;
      T omega1 ;
      T omega2 ;
      plint nx ;
      plint ny ;
      plint nz ;
      T G ;
      T force_fluid1 ;
      T force_fluid2 ;
      T Gads_f1_s1 ;
      T Gads_f1_s2 ;
	  std::string fNameIn ;
      plint maxIter ; //max_no_it
      plint saveIter ;
      plint vtkIter ;
      plint statIter;
      T convergence ;
      T tempHold ;
      plint startNum ;
    //  plint runnum;// no of density simulations to run
      T nx1f1 ; //fluid configuration
      T nx2f1 ;
      T ny1f1 ;
      T ny2f1 ;
      T nz1f1 ;
      T nz2f1 ;
      T nx1f2 ;
      T nx2f2 ;
      T ny1f2 ;
      T ny2f2 ;
      T nz1f2 ;
      T nz2f2 ;
	  T rho1  ;
	  T rho2  ;
      T rhoNoFluid ;
      T cap_factor ;
	  T rho_fluid2_min ;
	  T rho_fluid2_max ;
	  T rho_fluid2_step ;
	  T rho_fluid1_val ;
	  bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; //periodicity



  string xmlFname;
  try {
      global::argv(1).read(xmlFname);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: "
              << (std::string)global::argv(0) << " input-file.xml" << std::endl;
        return -1;
    }


    // 2. Read input parameters from the XML file.
      cout << "Reading inputs from xml file \n";
      try {
          XMLreader document(xmlFname);
	document["load"]["new_file"].read(new_file);
    document["geometry"]["file_geom"].read(fNameIn);
    document["geometry"]["size"]["x"].read(nx);
    document["geometry"]["size"]["y"].read(ny);
    document["geometry"]["size"]["z"].read(nz);

    document["init"]["fluid1"]["x1"].read(nx1f1);
    document["init"]["fluid1"]["x2"].read(nx2f1);
    document["init"]["fluid1"]["y1"].read(ny1f1);
    document["init"]["fluid1"]["y2"].read(ny2f1);
    document["init"]["fluid1"]["z1"].read(nz1f1);
    document["init"]["fluid1"]["z2"].read(nz2f1);

    document["init"]["fluid2"]["x1"].read(nx1f2);
    document["init"]["fluid2"]["x2"].read(nx2f2);
    document["init"]["fluid2"]["y1"].read(ny1f2);
    document["init"]["fluid2"]["y2"].read(ny2f2);
    document["init"]["fluid2"]["z1"].read(nz1f2);
    document["init"]["fluid2"]["z2"].read(nz2f2);

    document["fluids"]["Gc"].read(G);
    document["fluids"]["omega1"].read(omega1);
    document["fluids"]["omega2"].read(omega2);
    document["fluids"]["force1"].read(force_fluid1);
    document["fluids"]["force2"].read(force_fluid2);
    document["fluids"]["G_ads_f1_s1"].read(Gads_f1_s1);
    document["fluids"]["G_ads_f1_s2"].read(Gads_f1_s2);
	
    document["fluids"]["num_pc"].read(rho_fluid2_step);
    document["fluids"]["rho1_i"].read(rho_fluid1_val);
    document["fluids"]["rho2_i"].read(rho_fluid2_max);
 //   document["fluids"]["rho1_f"].read(rho1_f);
    document["fluids"]["rho2_f"].read(rho_fluid2_min);
	document["fluids"]["rho1"].read(rho1);
	document["fluids"]["rho2"].read(rho2);
    document["fluids"]["rho_d"].read(rhoNoFluid);
	
    document["output"]["conv"].read(convergence);
    document["output"]["out_f"].read(fNameOut);
    document["output"]["it_max"].read(maxIter);
    document["output"]["it_con"].read(statIter);
  // document["output"]["it_info"].read(statIter);
    document["output"]["it_gif"].read(saveIter);
    document["output"]["it_vtk"].read(vtkIter);

    document["geometry"]["per"]["fluid1"]["x"].read(px_f1);
    document["geometry"]["per"]["fluid1"]["y"].read(py_f1);
    document["geometry"]["per"]["fluid1"]["z"].read(pz_f1);
    document["geometry"]["per"]["fluid2"]["x"].read(px_f2);
    document["geometry"]["per"]["fluid2"]["y"].read(py_f2);
    document["geometry"]["per"]["fluid2"]["z"].read(pz_f2);

    }
      catch (PlbIOException& exception) {
          pcout << exception.what() << std::endl;
          pcout << exception.what() << std::endl;
          return -1;
      }
		  
	plint runnum = ((rho_fluid2_max - rho_fluid2_min) / rho_fluid2_step)+1;
    global::directories().setOutputDir(fNameOut);	
    T rho_fluid1[runnum] ;
    T rho_fluid2[runnum];
	//T mu1[runnum];
	//T mu2[runnum];
	T M[runnum];
	T Ca[runnum];
	// T Mom1[runnum];
	// T Mom2[runnum];
	// T Kr1[runnum];
	// T Kr2[runnum];
	// T Mom1_high;
	// T Mom2_high;
	T deltaP[runnum];
	// T kr1[runnum];
	// T kr2[runnum];
	T k1_high;
	T k2_high;
	T meanRho1;
	T meanRho2;
	T mu1;
	T mu2;
	T M1;
	T Ca1;
	T rho_F1;
	T rho_F2;
	// T Force1;
	// T Force2;
	T mean_U1[runnum];
	T mean_U2[runnum];
	T mean_rho1[runnum];
	T mean_rho2[runnum];
	// T diff_U1[runnum];
	// T diff_U2[runnum];
	// T diff_rho1[runnum];
	// T diff_rho2[runnum];
	
	const T delta_rho = 0.00005;
	T rho_low = 2 - delta_rho;	  
        	

		
	for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid2[readnum] = rho_fluid2_max - (readnum-1)* rho_fluid2_step;
		}
	for (plint readnum = 1; readnum <= runnum; ++readnum) {
        rho_fluid1[readnum] = rho_fluid1_val;
	}
	
    const T nu1 = ((T)1 / omega1 - 0.5) / DESCRIPTOR<T>::invCs2;
    const T nu2 = ((T)1 / omega2 - 0.5) / DESCRIPTOR<T>::invCs2;
	

    // Use regularized BGK dynamics to improve numerical stability (but note that BGK dynamics works well too).
    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega2));
    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1));

    lattice_fluid2.periodicity().toggle(0, px_f2);
    lattice_fluid1.periodicity().toggle(0, px_f1);
    lattice_fluid2.periodicity().toggle(1, py_f2);
    lattice_fluid1.periodicity().toggle(1, py_f1);
    lattice_fluid2.periodicity().toggle(2, pz_f2);
    lattice_fluid1.periodicity().toggle(2, pz_f1);

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
		deltaP[readnum]=(rho_fluid1[readnum]-rho_fluid2[readnum])/3;
		pcout << "Run number = " << readnum << endl;
        pcout << "Rho_no_1 = " << rho_fluid1[readnum] << endl;
        pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
    }

    pcout << "Cap Multiplier = " << cap_factor << endl;
    pcout << "Reading the geometry file." << endl;

     MultiScalarField3D<int> geometry(nx, ny, nz);
	 readGeometry(fNameIn, fNameOut, geometry);
	 
    util::ValueTracer<T> converge1(1.0, lattice_fluid2.getNx(), convergence); // convergence parameters velocity/size/threshold
    util::ValueTracer<T> converge2(1.0, lattice_fluid2.getNx(), convergence);

	
    for (plint runs = 1; runs <= runnum; ++runs) { // Loop simulations with varying saturation
	
	 pcout << "Run number = " << runs << endl;

		
        if (runs > 1) {
            pcout << "Using previous Geometry  " << endl;
        }
        else {
			
            PorousMediaSetup(lattice_fluid1, lattice_fluid2, geometry, rhoNoFluid, rho1, rho2, Gads_f1_s1, Gads_f1_s2, force_fluid1, force_fluid2, nx1f1, nx2f1, ny1f1, ny2f1, nz1f1, nz2f1, nx1f2, nx2f2, ny1f2, ny2f2, nz1f2, nz2f2, runs, new_file);
			
	   }
		


        T meanJ_2 = 10000.0;
        T meanJ_old2;
        T meanJ_1 = 10000.0;
        T meanJ_old1;
        T convergemeanJ1;
        T convergemeanJ2;
        T convergedmeanJ1;
        T convergedmeanJ2;
		

		
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


			converge1.takeValue(getStoredAverageDensity(lattice_fluid1), true); //check for convergence
			converge2.takeValue(getStoredAverageDensity(lattice_fluid2), true); //check for convergence
            if (iT % statIter == 0) {
				mean_rho1[runs] = getStoredAverageDensity<T>(lattice_fluid1) ;
				mean_rho2[runs] = getStoredAverageDensity<T>(lattice_fluid2);
                pcout << "Iteration:  " << iT << endl;
                pcout << "Average density fluid one = " << mean_rho1[runs] << endl;
                pcout << "Average density fluid two = " << mean_rho2[runs] << endl << endl;


                if ((converge1.hasConverged()) && (converge2.hasConverged())) {
					pcout << "Converge 1 and Converge 2 converged for iteration:" << iT << endl; 
                    checkconv = 1;
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
					writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
					
				// Calculate velocity here for both fluids in x-direction and pcout
					
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				mean_U1[runs] = meanU1;
				mean_U2[runs] = meanU2;

				
				T rho_F1=rho_fluid1[runs];
				T rho_F2=rho_fluid2[runs];
				
				// calculate capillary number & dynamic viscosity ratio
				computeRatios(lattice_fluid1, lattice_fluid2, rho_F1, rho_F2, nu1, nu2, Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1, meanRho1, meanRho2, M1, Ca1);
				M[runs]=M1;
				Ca[runs]=Ca1;
						
                }


                if ((convergemeanJ1 < convergedmeanJ1) && (convergemeanJ2 < convergedmeanJ2)) {
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    pcout << "Simulation has converged" << endl;
                    writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                    writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
					
				// Calculate velocity here for both fluids in x-direction and pcout
					
				
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				mean_U1[runs] = meanU1;
				mean_U2[runs] = meanU2;
				
		
				// calculate capillary number & dynamic viscosity ratio
				computeRatios(lattice_fluid1, lattice_fluid2, rho_F1, rho_F2, nu1, nu2, Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1, meanRho1, meanRho2, M1, Ca1);
				M[runs]=M1;
				Ca[runs]=Ca1;
				
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
				
				T meanU1 = computeVelocity1(lattice_fluid1, nu1);
				T meanU2 = computeVelocity2(lattice_fluid2, nu2);
				mean_U1[runs] = meanU1;
				mean_U2[runs] = meanU2;
				
				
				T rho_F1=rho_fluid1[runs];
				T rho_F2=rho_fluid2[runs];
				
				// calculate capillary number & dynamic viscosity ratio
				computeRatios(lattice_fluid1, lattice_fluid2, rho_F1, rho_F2, nu1, nu2, Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1, meanRho1, meanRho2, M1, Ca1);
				M[runs]=M1;
				Ca[runs]=Ca1;
				
				
            }
			
            initializeAtEquilibrium(lattice_fluid1, Box3D(1, 2, 0, ny-1, 0, nz-1), rho_fluid1[runs], Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid2, Box3D(1, 2, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid1, Box3D(nx - 2, nx-1, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
            initializeAtEquilibrium(lattice_fluid2, Box3D(nx - 2, nx-1, 0, ny-1, 0, nz-1), rho_fluid2[runs], Array<T, 3>(0., 0., 0.));
			 
            lattice_fluid1.initialize();
            lattice_fluid2.initialize();		
        }
    }
	//  outputting variables
	
	saveBinaryBlock(lattice_fluid1, "lattice_fluid1.dat");
	saveBinaryBlock(lattice_fluid2, "lattice_fluid2.dat");
	
	
	
	std::string outDir = fNameOut + "/";
	std::string output = outDir + fNameIn + "_output.dat";
	t = clock() - t;
	pcout << "Simulation took seconds:" << ((float)t)/CLOCKS_PER_SEC << std::endl;
	plb_ofstream ofile(output.c_str());
	ofile << "Outputs" << "\n\n";
	ofile << "Simulation took seconds:" << ((float)t)/CLOCKS_PER_SEC << endl;
	

	
	for (plint runs = 1; runs <= runnum; ++runs) {
	
	
	pcout << "Run    = " << runs            << std::endl;
	pcout << "Pressure difference =  " << deltaP[runs] << std::endl;
	pcout << "Viscosity ratio =  " << M[runs] << std::endl;
	pcout << "Capillary number =  " << Ca[runs] << std::endl;
	
	
	ofile << "Run = " << runs << "\n" << endl;
	ofile << "Pressure difference = " << deltaP[runs] <<"\n" << endl;
	ofile << "Viscosity ratio =  " << M[runs] <<"\n" << endl;
	ofile << "Capillary number =  " << Ca[runs] <<"\n" << endl;
	
	}
	
}
