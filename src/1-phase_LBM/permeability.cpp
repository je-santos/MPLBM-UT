/* Modified from file by Wim Degruyter */

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


using namespace plb;

typedef double T;
//#define DESCRIPTOR descriptors::D3Q19Descriptor
#define DESCRIPTOR descriptors::MRTD3Q19Descriptor

// This function object returns a zero velocity, and a pressure which decreases
//   linearly in x-direction. It is used to initialize the particle populations.
class PressureGradient {
public:
  PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
  { }
  void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
  {
    velocity.resetToZero();
    density = (T)1 - deltaP*DESCRIPTOR<T>::invCs2 / (T)(nx-1) * (T)iX;
  }
private:
  T deltaP;
  plint nx;
};

// This function grabs the appropiate geometry for the single-phase simulation
void readGeometry(std::string fNameIn, std::string fNameOut,
  MultiScalarField3D<int>& geometry, plint run, plint runnum, bool vtk_out,
  std::string GeometryName)
{
  const plint nx = geometry.getNx();
  const plint ny = geometry.getNy();
  const plint nz = geometry.getNz();
  plint run_diff = ((runnum - 1)/2)+1;
  std::string fNameIn_temp1 = fNameIn + GeometryName;
  std::string fNameIn_temp = "0";

  pcout  << "\n" << " Run: "<< run << std::endl;

  Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);

  // Selects between the original geometry or the fluid 1,2 final config from multiphase sim

  if (run == 1) {  // original geometry - absolute permeability
    fNameIn_temp = fNameIn_temp1 + ".dat";
    pcout  << "Running absolute permeability "<< std::endl;
  }

  if (run > run_diff) { // Fluid 1 : Krnw
    const plint runner = run - run_diff;
    fNameIn_temp =  fNameIn + "f1_for_kr_" + patch::to_string(runner)+".dat";
    pcout  << "Running kr_f1 "<< std::endl;
  }

  if ( run > 1 && run < (run_diff+1) ) { // Fluid 2 : Krw
    fNameIn_temp =  fNameIn + "f2_for_kr_"+ patch::to_string(run-1)+".dat";
    pcout  << "Running kr_f2 "<< std::endl;
  }


  pcout  << "The geometry name is  "<< fNameIn_temp << std::endl;

  std::unique_ptr<plb::MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
  plb_ifstream geometryFile(fNameIn_temp.c_str());

  for (plint iX=3; iX<nx-4; ++iX) {
    if (!geometryFile.is_open()) {
      pcout << "Error: could not open the geometry file " << fNameIn_temp << std::endl;
      exit(EXIT_FAILURE);
    }

    geometryFile >> *slice;
    copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
  }

  if  (vtk_out == true) {
    VtkImageOutput3D<T> vtkOut(createFileName("PorousMedium", run, 6), 1.0);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);
  }


 // code to create .st file. Uncomment if needed
  //{
    //std::auto_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int,T>(geometry, geometry.getBoundingBox());
    //std::vector<T> isoLevels;
    //isoLevels.push_back(0.5);
    //typedef TriangleSet<T>::Triangle Triangle;
    //std::vector<Triangle> triangles;
    //Box3D domain = floatTags->getBoundingBox().enlarge(-1);
    //domain.x0++;
    //domain.x1--;
    //isoSurfaceMarchingCube(triangles, *floatTags, isoLevels, domain);
    //TriangleSet<T> set(triangles);
    //std::string outDir = fNameOut + "/";
    //set.writeBinarySTL(outDir + "porousMedium.stl");
  //}
}

void porousMediaSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
  MultiScalarField3D<int>& geometry, T deltaP)
  {
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    pcout << "Definition of inlet/outlet." << std::endl;

    Box3D inlet (0,0, 1,ny-2, 1,nz-2);
    boundaryCondition->addPressureBoundary0N(inlet, lattice);
    setBoundaryDensity(lattice, inlet, (T) 1.);

    Box3D outlet(nx-1,nx-1, 1,ny-2, 1,nz-2);
    boundaryCondition->addPressureBoundary0P(outlet, lattice);
    setBoundaryDensity(lattice, outlet, (T) 1. - deltaP*DESCRIPTOR<T>::invCs2);

    // Where "geometry" evaluates to 1, use bounce-back.
    defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
    // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
    defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

    //   pcout << "Initialization of rho and u." << std::endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(),
                                PressureGradient(deltaP, nx) );

    lattice.initialize();
    delete boundaryCondition;
  }

  void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, plint run)
  {
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo");

    // Write velocity-norm at x=1.
    imageWriter.writeScaledGif(createFileName("ux_inlet", run, 6),
    *computeVelocityNorm(lattice, Box3D(2,2, 0,ny-1, 0,nz-1)),
    imSize, imSize );

    // Write velocity-norm at x=nx/2.
    imageWriter.writeScaledGif(createFileName("ux_half", run, 6),
    *computeVelocityNorm(lattice, Box3D(nx/2,nx/2, 0,ny-1, 0,nz-1)),
    imSize, imSize );
  }

  void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, plint run)
  {
    VtkImageOutput3D<T> vtkOut(createFileName("vtk_vel", run, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
  }

  void computePermeability(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T nu, T deltaP, Box3D domain, T& perm, T& meanU)
  {

    // Compute only the x-direction of the velocity (direction of the flow).
    plint xComponent = 0;
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();
    Box3D domain1(0, nx-1, 0, ny-1, 0, nz-1);

    meanU = computeAverage(*computeVelocityComponent(lattice, domain1, xComponent));

    pcout << "Average velocity     = " << meanU                    << std::endl;
    //pcout << "Lattice viscosity nu = " << nu                       << std::endl;
    //pcout << "Grad P               = " << deltaP/(T)(nx-1)         << std::endl;
    perm = nu*meanU / (deltaP/(T)(nx-1));
    pcout << "Permeability         = " << perm 						<< std::endl;
    //  return meanU;
  }

  int main(int argc, char **argv)
  {
    plbInit(&argc, &argv);


    std::string fNameIn ;
    std::string fNameOut;

    plint nx;
    plint ny;
    plint nz;
    T deltaP ;
    T Run;
	  bool nx_p, ny_p, nz_p;
	  bool vtk_out;
    std::string GeometryName ;
    plint maxT;
	  T conv;

    std::string xmlFname;
    try {
      global::argv(1).read(xmlFname);
    }
    catch (PlbIOException& exception) {
      pcout << "Wrong parameters; the syntax is: "
      << (std::string) global::argv(0) << " input-file.xml" << std::endl;
      return -1;
    }

    // 2. Read input parameters from the XML file.
    pcout << "Reading inputs from xml file \n";
    try {
      XMLreader document(xmlFname);
      document["geometry"]["file_geom"].read(GeometryName);

      document["geometry"]["size"]["x"].read(nx);
      document["geometry"]["size"]["y"].read(ny);
      document["geometry"]["size"]["z"].read(nz);
	    document["geometry"]["per"]["x"].read(nx_p);
      document["geometry"]["per"]["y"].read(ny_p);
      document["geometry"]["per"]["z"].read(nz_p);


      document["folder"]["out_f"].read(fNameOut);
      document["folder"]["in_f"].read(fNameIn);

      document["simulations"]["press"].read(deltaP);
      document["simulations"]["num"].read(Run);
      document["simulations"]["iter"].read(maxT);
	    document["simulations"]["conv"].read(conv);
	    document["simulations"]["vtk_out"].read(vtk_out);

    }
    catch (PlbIOException& exception) {
      pcout << exception.what() << std::endl;
      pcout << exception.what() << std::endl;
      return -1;
    }


    std::string inputF= fNameIn;
    global::directories().setOutputDir(fNameOut+"/");
    global::directories().setInputDir(inputF+"/");

    const T omega = 1.0;
    const T nu    = ((T)1/omega- (T)0.5)/DESCRIPTOR<T>::invCs2;
    const plint runnum = Run;
    plint run_diff = ((runnum - 1)/2)+1;

    T perm[runnum];
    T meanU[runnum];
    T rel_perm[runnum];
    T Perm;
    T Vel;
    pcout << "Total simulations: " << runnum << std::endl;
    pcout << "The convergence threshold is: " << conv << " %" << std::endl;


    for (plint run = 1; run <= runnum; ++run) {

      MultiBlockLattice3D<T,DESCRIPTOR> lattice( nx,ny,nz,
                                        new BGKdynamics<T,DESCRIPTOR>(omega) );
      // Switch off periodicity.
      //lattice.periodicity().toggleAll(false);

      lattice.periodicity().toggle(0, nx_p);
      lattice.periodicity().toggle(1, ny_p);
      lattice.periodicity().toggle(2, nz_p);

      MultiScalarField3D<int> geometry(nx,ny,nz);
      readGeometry(fNameIn, fNameOut, geometry, run, runnum, vtk_out, GeometryName);

      porousMediaSetup(lattice, createLocalBoundaryCondition3D<T,DESCRIPTOR>(),
                       geometry, deltaP);


      pcout << "Simulation begins" << std::endl;
      plint iT=0;
      T new_avg_f, old_avg_f, relE_f1;
      lattice.toggleInternalStatistics(false);

      for (;iT<maxT; ++iT) {


        if (iT % 250 == 0 && iT > 0) {

          lattice.toggleInternalStatistics(true);
          pcout << "Iteration " << iT   << std::endl;
          pcout << "-----------------"  << std::endl;
          lattice.collideAndStream();
          new_avg_f = getStoredAverageEnergy(lattice);
          lattice.toggleInternalStatistics(false);
          relE_f1 = std::fabs(old_avg_f-new_avg_f)*100/old_avg_f;
          pcout << "Relative difference of Energy: " << setprecision(3)
          << relE_f1 <<" %"<<std::endl;
          pcout << "The preliminary permeability is: " <<std::endl;
          computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(), Perm, Vel);
          pcout << "**********************************************" <<std::endl;
          if ( relE_f1<conv ){
            break;
          }
          old_avg_f = new_avg_f; // store new properties
        }
      }

      pcout << "End of simulation at iteration " << iT << " for Run "<< run << std::endl;

      //   pcout << "Permeability:" << std::endl;
      computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(), Perm, Vel);

	    writeGifs(lattice,iT,run);
	    std::string outDir = fNameOut + "/";
      std::string vel_name = outDir + GeometryName + "_vel.dat";
      plb_ofstream ofile3( vel_name.c_str() );
      ofile3 << setprecision(1) <<*computeVelocity(lattice) << endl;
        
      //std::string rho_name = outDir + GeometryName + "_rho.dat";
      //plb_ofstream ofile4( rho_name.c_str() );
      //ofile4 << setprecision(10) <<*computeDensity(lattice) << endl;

      perm[run]=Perm;
      meanU[run]=Vel;

      rel_perm[run]=perm[run]/perm[1];
      if (run == 1) {
        pcout << "Absolute Permeability   = " << perm[run]         << std::endl;
      }
      pcout << "Relative Permeability = " << rel_perm[run] << std::endl;

	 if  (vtk_out == true) {
      pcout << "Writing VTK file ..." << std::endl;
      writeVTK(lattice, iT, run);
	 }
    }

    pcout << "Printing outputs" << std::endl;
    std::string outDir = fNameOut + "/";
    std::string output = outDir + "relPerm&vels.txt";
    plb_ofstream ofile(output.c_str());
    ofile << "Outputs" << "\n\n";
    ofile << "Krw from run: " << "\n" << "Krnw from run: " << (run_diff+1) << std::endl;
    for (plint runs = 1; runs <= runnum; ++runs) {


      ofile << "Run   = " << runs        << std::endl;
      if (runs == 1) {
        ofile << "Absolute Permeability   = " << perm[runs]         << std::endl;

      }

      ofile << "Relative Permeability   = " << rel_perm[runs]         << std::endl;
      ofile << "Mean Velocity   = " << meanU[runs]         << std::endl;

    }
  }
