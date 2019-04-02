/*************************
Palabos Shan-Chen Multiphase Simulator
Capillary Pressure Increments
author: Javier E. Santos
Based on the code by Dr. Christopher Landry
University of Texas at Austin


*************************/

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

typedef double T;
#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor // Use a grid which additionally to the f's stores two variables for the external force term.
static std::string outputDir("./tmp/");

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.

struct Param
{
  T omega1, omega2;
  T nu1, nu2;                               // Kinematic viscosity.
  //T lx, ly, lz;                       // Size of computational domain, in physical units.
  plint it_max, it_con, it_info, it_gif, it_vtk;
  //plint resolution;                   // Number of lattice nodes along a reference length.
  plint nx, ny, nz;                   // Grid resolution of bounding box.
  plint nx1f1,nx2f1,ny1f1,ny2f1,nz1f1,nz2f1;
  plint nx1f2,nx2f2,ny1f2,ny2f2,nz1f2,nz2f2;
  plint num_pc;
  T rho1_i, rho2_i, rho1_f, rho2_f, rho1, rho2, rho_d;
  T Gc;
  T force_fluid1, force_fluid2;
  T Gads_f1_s1, Gads_f1_s2;
  bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; //periodicity
  T *rho_fluid1;
  T *rho_fluid2;
  std::string geometry_fname,output_folder;

  Param()
  { }

  Param(std::string xmlFname)
  {
    XMLreader document(xmlFname);
    document["geometry"]["file_geom"].read(geometry_fname);
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

    document["fluids"]["Gc"].read(Gc);
    document["fluids"]["omega1"].read(omega1);
    document["fluids"]["omega2"].read(omega2);
    document["fluids"]["force1"].read(force_fluid1);
    document["fluids"]["force2"].read(force_fluid2);
    document["fluids"]["G_ads_f1_s1"].read(Gads_f1_s1);
    document["fluids"]["G_ads_f1_s2"].read(Gads_f1_s2);
    document["fluids"]["num_pc"].read(num_pc);
    document["fluids"]["rho1_i"].read(rho1_i);
    document["fluids"]["rho2_i"].read(rho2_i);
    document["fluids"]["rho1_f"].read(rho1_f);
    document["fluids"]["rho2_f"].read(rho2_f);
    document["fluids"]["rho1"].read(rho1);
    document["fluids"]["rho2"].read(rho2);
    document["fluids"]["rho_d"].read(rho_d);

    document["output"]["out_f"].read(output_folder);
    document["output"]["it_max"].read(it_max);
    document["output"]["it_con"].read(it_con);
    document["output"]["it_info"].read(it_info);
    document["output"]["it_gif"].read(it_gif);
    document["output"]["it_vtk"].read(it_vtk);

    document["geometry"]["per"]["fluid1"]["x"].read(px_f1);
    document["geometry"]["per"]["fluid1"]["y"].read(py_f1);
    document["geometry"]["per"]["fluid1"]["z"].read(pz_f1);
    document["geometry"]["per"]["fluid2"]["x"].read(px_f2);
    document["geometry"]["per"]["fluid2"]["y"].read(py_f2);
    document["geometry"]["per"]["fluid2"]["z"].read(pz_f2);


    pcout << "LBM constants are computed here\n";
    rho_fluid1= new T[num_pc];
    rho_fluid2= new T[num_pc];
    rho_fluid1[0]=rho1_i;
    rho_fluid1[num_pc-1]=rho1_f;
    for (plint i=1; i<num_pc; i++){
      rho_fluid1[i]=rho1_i+(rho1_f-rho1_i)/(num_pc-1)*i;
    }
    const T nu1 = ((T)1 / omega1 - 0.5) / DESCRIPTOR<T>::invCs2;
    const T nu2 = ((T)1 / omega2 - 0.5) / DESCRIPTOR<T>::invCs2;
  } // end of reading function
}; // end of struct

Param param;

void print_variables()
{
  pcout << "Convergence = " << param.it_con << endl;
  pcout << "nx = " << param.nx << endl;
  pcout << "ny = " << param.ny << endl;
  pcout << "nz = " << param.nz << endl;
  pcout << "Gc = " << param.Gc << endl;
  pcout << "force fluid 1 = " << param.force_fluid1 << endl;
  pcout << "force fluid 2 = " << param.force_fluid2 << endl;
  pcout << "G_ads_1 fluid 1 = " << param.Gads_f1_s1 << endl;
  pcout << "G_ads_1 fluid 2 = " << -param.Gads_f1_s1 << endl;
  pcout << "G_ads_2 fluid 1 = " << param.Gads_f1_s2 << endl;
  pcout << "G_ads_2 fluid 2 = " << -param.Gads_f1_s2 << endl;
  pcout << "nz2f2  = " << param.nz2f2 << endl;


  for (plint readnum = 1; readnum <= param.num_pc; ++readnum)
    {
    pcout << "Rho_no_1 = " << param.rho_fluid1[readnum] << endl;
    }
}

int porousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, MultiScalarField3D<int>& geometry){

      pcout << "Definition of the geometry." << std::endl;

      // Where "geometry" evaluates to 1, use bounce-back (wettability 1).
    defineDynamics(lattice_fluid1, geometry, new BounceBack<T,DESCRIPTOR>( param.Gads_f1_s1), 1);
    defineDynamics(lattice_fluid1, geometry, new BounceBack<T,DESCRIPTOR>(-param.Gads_f1_s1), 1);
    // Where "geometry" evaluates to 3, use bounce-back (wettability 2).
    defineDynamics(lattice_fluid1, geometry, new BounceBack<T,DESCRIPTOR>( param.Gads_f1_s2), 2);
    defineDynamics(lattice_fluid1, geometry, new BounceBack<T,DESCRIPTOR>(-param.Gads_f1_s2), 2);
    // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
    defineDynamics(lattice_fluid1, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);
    defineDynamics(lattice_fluid1, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

    Array<T, 3> zeroVelocity(0., 0., 0.);
    // Initialize  uniform density for target saturation
    pcout << "Initializing Fluids" << endl;

  //initializeAtEquilibrium(lattice_fluid2, Box3D(param.nx1f2, param.nx2f2, param.ny1f2, param.ny2f2, param.nz1f2, param.nz2f2), param.rho_fluid2, zeroVelocity);
  //initializeAtEquilibrium(lattice_fluid1, Box3D(param.nx1f2, param.nx2f2, param.ny1f2, param.ny2f2, param.nz1f2, param.nz2f2), param.rho_d, zeroVelocity);
  //initializeAtEquilibrium(lattice_fluid1, Box3D(param.nx1f1, param.nx2f1, param.ny1f1, param.ny2f1, param.nz1f1, param.nz2f1), param.rho_fluid1, zeroVelocity);
  //initializeAtEquilibrium(lattice_fluid2, Box3D(param.nx1f1, param.nx2f1, param.ny1f1, param.ny2f1, param.nz1f1, param.nz2f1), param.rho_d, zeroVelocity);

    setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., param.force_fluid1, 0.));
    setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., param.force_fluid2, 0.));

    lattice_fluid1.initialize();
    lattice_fluid2.initialize();

      //for (plint itX = 0; itX < param.nx; ++itX) {
      //  for (plint itY = 0; itY < param.ny; ++itY) {
      //      for (plint itZ = 0; itZ < param.nz; ++itZ) {
      //          if (geometry.get(itX, itY, itZ) == 1) {
      //              defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>( param.Gads_f1_s1));
      //              defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-param.Gads_f1_s1));
      //          }
      //          if (geometry.get(itX, itY, itZ) == 2) {
      //              defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
      //              defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
      //          }
      //          if (geometry.get(itX, itY, itZ) == 3) {
      //              defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>( param.Gads_f1_s2));
      //              defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-param.Gads_f1_s2));
      //          }

              //  }
          //  }
      //  }
  //  }

  return 0;
  }


int create_geometry()
{



  //vector<MultiBlockLattice3D<T, DESCRIPTOR>*> blockLattices;
  //blockLattices.push_back(&lattice_fluid2);
  //blockLattices.push_back(&lattice_fluid1);
  //std::vector<T> constOmegaValues;
  //constOmegaValues.push_back(param.omega2);
  //constOmegaValues.push_back(param.omega1);
  //plint processorLevel = 1;
  //integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D<T,
//    DESCRIPTOR>(param.Gc, constOmegaValues),

//  Box3D(0, param.nx-1, 0, param.ny-1, 0, param.nz-1), blockLattices, processorLevel);
//  MultiScalarField3D<int> geometry(param.nx, param.ny, param.nz);
//  plb_ifstream geometryFile(param.geometry_fname.c_str());

//  if (!geometryFile.is_open()) {
//      pcout << "Error: could not open geometry file " << param.geometry_fname << endl;
//      return -1;
//  }
//  geometryFile >> geometry;

return 0;
}



int main(int argc, char* argv[])
{

  plbInit(&argc, &argv);
  global::directories().setOutputDir(outputDir);
  // The try-catch blocks catch exceptions in case an error occurs,
  // and terminate the program properly with a nice error message.
  // 1. Read command-line parameter: the input file name.
  string xmlFileName;
  try {
      global::argv(1).read(xmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: "
              << (std::string)global::argv(0) << " input-file.xml" << std::endl;
        return -1;
    }
    // 2. Read input parameters from the XML file.
      cout << "Reading inputs from xml file \n";
      try {
          param = Param(xmlFileName);
      }
      catch (PlbIOException& exception) {
          pcout << exception.what() << std::endl;
          return -1;
      }


  cout << "Print variables function: \n";
  print_variables();

  cout << "Let's read the geometry: \n";

  pcout << "Creation of the lattice." << endl;
  MultiBlockLattice3D<T,DESCRIPTOR> lattice_fluid1(param.nx,param.ny,param.nz,
          new ExternalMomentRegularizedBGKdynamics<T,DESCRIPTOR>(param.omega1));
  MultiBlockLattice3D<T,DESCRIPTOR> lattice_fluid2(param.nx,param.ny,param.nz,
          new ExternalMomentRegularizedBGKdynamics<T,DESCRIPTOR>(param.omega2));
  // Switch off periodicity.
  lattice_fluid1.periodicity().toggleAll(false);
  lattice_fluid2.periodicity().toggleAll(false);

  pcout << "Creating geometry file\n";

  MultiScalarField3D<int> geometry(param.nx, param.ny, param.nz);
  Box3D sliceBox(0,0, 0,param.ny-1, 0,param.nz-1);
  std::auto_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
  plb_ifstream geometryFile(param.geometry_fname.c_str());
  for (plint iX=0; iX<param.nx-1; ++iX) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << param.geometry_fname << std::endl;
            exit(EXIT_FAILURE);
          }
      geometryFile >> *slice;
      copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,param.ny-1, 0,param.nz-1));
    }

    VtkImageOutput3D<T> vtkOut("porousMedium", 1.0);
    vtkOut.writeData<float>(*copyConvert<int,T>(geometry, geometry.getBoundingBox()), "tag", 1.0);

    cout << "Let's set up the PM: \n";
     porousMediaSetup(lattice_fluid1, lattice_fluid2, geometry);
cout << "Running the code \n";
plint iT = 0;
for (iT=0; iT<=100; ++iT){

     lattice_fluid1.collideAndStream(); // check order
    lattice_fluid2.collideAndStream();
}
    //initializeAtEquilibrium(lattice_fluid1, Box3D(1, 2, 0, ny, 0, nz), rho_fluid1[runs] * cap_factor, Array<T, 3>(0., 0., 0.));
    //        initializeAtEquilibrium(lattice_fluid2, Box3D(1, 2, 0, ny, 0, nz), rhoNoFluid, Array<T, 3>(0., 0., 0.));
    //        initializeAtEquilibrium(lattice_fluid1, Box3D(nx - 2, nx-1, 0, ny, 0, nz), rhoNoFluid, Array<T, 3>(0., 0., 0.));
    //        initializeAtEquilibrium(lattice_fluid2, Box3D(nx - 2, nx-1, 0, ny, 0, nz), rho_fluid2[runs], Array<T, 3>(0., 0., 0.));

    //        lattice_fluid1.initialize();
    //        lattice_fluid2.initialize();

return 0;
}
