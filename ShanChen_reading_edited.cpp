/*************************
Palabos Shan-Chen Multiphase Simulator
Capillary Pressure Increments
author: Javier E. Santos & Abhishek Bihani
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
    T nu1,nu2;                               // Kinematic viscosity.
    //T lx, ly, lz;                       // Size of computational domain, in physical units.
    plint it_max, it_con, it_info, it_gif, it_vtk;
    //plint resolution;                   // Number of lattice nodes along a reference length.
    plint nx, ny, nz;                   // Grid resolution of bounding box.
    plint nx1f1,nx2f1,ny1f1,ny2f1,nz1f1,nz2f1;
    plint nx1f2,nx2f2,ny1f2,ny2f2,nz1f2,nz2f2;
    plint num_pc;
    T rho1_i, rho2_i, rho1_f, rho2_f, rho1, rho2, rhoNoFluid;
    T G;
    T force_fluid1, force_fluid2;
    T Gads_f1_s1, Gads_f1_s2;
    bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; //periodicity
    T *rho_fluid1;
    T *rho_fluid2;
    T nu1, nu2;
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

        document["fluids"]["omega1"].read(omega1);
        document["fluids"]["omega2"].read(omega2);
        document["fluids"]["force1"].read(force_fluid1);
        document["fluids"]["force2"].read(force_fluid2);
        document["fluids"]["G_ads_f1_s1"].read(Gads_f1_s1);
        document["fluids"]["G_ads_f1_s2"].read(Gads_f1_s2);
        document["fluids"]["Gc"].read(G);
        document["fluids"]["num_pc"].read(num_pc);
        document["fluids"]["rho1_i"].read(rho1_i);
        document["fluids"]["rho2_i"].read(rho2_i);
        document["fluids"]["rho1_f"].read(rho1_f);
        document["fluids"]["rho2_f"].read(rho2_f);
        document["fluids"]["rho1"].read(rho1);
        document["fluids"]["rho2"].read(rho2);
        document["fluids"]["rho_d"].read(rhoNoFluid);

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
//}

  //  void computeLBparameters()
//    {
    pcout << "LBM constants are computed here/n";
    rho_fluid1= new T[num_pc]
    rho_fluid2= new T[num_pc]
    rho_fluid1[0]=rho1_i;
    rho_fluid1[num_pc-1]=rho1_f;
    rho_fluid2[0]=rho2_i;
    rho_fluid2[num_pc-1]=rho2_f;
    for (plint i=1; i<num_pc; i++){
      rho_fluid1[i]=rho1_i+(rho1_f-rho1_i)/(num_pc-1)*i;
      rho_fluid2[i]=rho2_i+(rho2_f-rho2_i)/(num_pc-1)*i;
    }
    const T nu1 = ((T)1 / omega1 - 0.5) / DESCRIPTOR<T>::invCs2;
    const T nu2 = ((T)1 / omega2 - 0.5) / DESCRIPTOR<T>::invCs2;
    }

  };

  Param param;
//  void dummy_func()
//  {
//    pcout << "nx" << param.nx;
//  }
void PorousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, MultiScalarField3D<int>& geometry)
{
//  plint nx = lattice_fluid2.getNx();
//  plint ny = lattice_fluid2.getNy();
//  plint nz = lattice_fluid2.getNz();
  pcout << "rho_fluid1[0]:" << param.rho_fluid1[0] << endl;
  pcout << "Definition of the geometry." << endl;
  for (plint itX = 0; itX < param.nx; ++itX) {
      for (plint itY = 0; itY < param.ny; ++itY) {
          for (plint itZ = 0; itZ < param.nz; ++itZ) {
              if (geometry.get(itX, itY, itZ) == 1) {
                  defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(param.Gads_f1_s1));
                  defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-param.Gads_f1_s1));
              }
              if (geometry.get(itX, itY, itZ) == 2) {
                  defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
                  defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
              }
              if (geometry.get(itX, itY, itZ) == 3) {
                  defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(param.Gads_f1_s2));
                  defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-param.Gads_f1_s2));
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
  VtkImageOutput3D<int> vtkOut(createFileName("Geometry", 1, 1), 1.);
  //vtkOut.writeData<int>(geometry, "Dynamics", 1.);
  vtkOut.writeData<int>(geometry, "Domain", 1.);

  Array<T, 3> zeroVelocity(0., 0., 0.);

  // Initialize  uniform density for target saturation
  pcout << "Initializing Fluids" << endl;

  initializeAtEquilibrium(lattice_fluid2, Box3D(param.nx1f2, param.nx2f2, param.ny1f2, param.ny2f2, param.nz1f2, param.nz2f2), param.rho_fluid2, zeroVelocity);
  initializeAtEquilibrium(lattice_fluid1, Box3D(param.nx1f2, param.nx2f2, param.ny1f2, param.ny2f2, param.nz1f2, param.nz2f2), param.rhoNoFluid, zeroVelocity);
  initializeAtEquilibrium(lattice_fluid1, Box3D(param.nx1f1, param.nx2f1, param.ny1f1, param.ny2f1, param.nz1f1, param.nz2f1), param.rho_fluid1, zeroVelocity);
  initializeAtEquilibrium(lattice_fluid2, Box3D(param.nx1f1, param.nx2f1, param.ny1f1, param.ny2f1, param.nz1f1, param.nz2f1), param.rhoNoFluid, zeroVelocity);

  //initializeAtEquilibrium(lattice_fluid1, Box3D(0, 4, 0, ny, 0, nz), rho_fluid1, zeroVelocity);
  //initializeAtEquilibrium(lattice_fluid2, Box3D(0, 4, 0, ny, 0, nz), rhoNoFluid, zeroVelocity);

  //initializeAtEquilibrium(lattice_fluid2, Box3D(5, nx, 0, ny, 0, nz), rho_fluid2, zeroVelocity);
  //initializeAtEquilibrium(lattice_fluid1, Box3D(5, nx, 0, ny, 0, nz), rhoNoFluid, zeroVelocity);

  setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., param.force_fluid1, 0.));
  setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., param.force_fluid2, 0.));

  lattice_fluid1.initialize();
  lattice_fluid2.initialize();

}

 void print_variables()
 {
   pcout << "Convergence = " << param.it_con << endl;
   pcout << "nx = " << param.nx << endl;
   pcout << "ny = " << param.ny << endl;
   pcout << "nz = " << param.nz << endl;
   pcout << "Gc = " << param.G << endl;
   pcout << "force fluid 1 = " << param.force_fluid1 << endl;
   pcout << "force fluid 2 = " << param.force_fluid2 << endl;
   pcout << "G_ads_1 fluid 1 = " << param.Gads_f1_s1 << endl;
   pcout << "G_ads_1 fluid 2 = " << -param.Gads_f1_s1 << endl;
   pcout << "G_ads_2 fluid 1 = " << param.Gads_f1_s2 << endl;
   pcout << "G_ads_2 fluid 2 = " << -param.Gads_f1_s2 << endl;
   pcout << "nz2f2  = " << param.nz2f2 << endl;
   pcout << "Rho_no_fluid = " << param.rhoNoFluid << endl;

   for (plint readnum = 1; readnum <= param.num_pc; ++readnum) {
       pcout << "Rho_no_1 = " << param.rho_fluid1[readnum] << endl;
       pcout << "Rho_no_2 = " << param.rho_fluid2[readnum] << endl;
   }

 }

void runProgram()
{
//  dummy_func();
  //computeLBparameters();
  cout << "test px_f1" << param.px_f1 << "\n";
  cout << "Read geometry: porousMediaSetup\n";
//  PorousMediaSetup();

  // Use regularized BGK dynamics to improve numerical stability (but note that BGK dynamics works well too).
  MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2(param.nx, param.ny, param.nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(param.omega2));
  MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1(param.nx, param.ny, param.nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(param.omega1));

  cout << "Set periodicity, omegas, etc\n";

  lattice_fluid2.periodicity().toggle(0, param.px_f2);
  lattice_fluid1.periodicity().toggle(0, param.px_f1);
  lattice_fluid2.periodicity().toggle(1, param.py_f2);
  lattice_fluid1.periodicity().toggle(1, param.py_f1);
  lattice_fluid2.periodicity().toggle(2, param.pz_f2);
  lattice_fluid1.periodicity().toggle(2, param.pz_f1);

  vector<MultiBlockLattice3D<T, DESCRIPTOR>*> blockLattices;
  blockLattices.push_back(&lattice_fluid2);
  blockLattices.push_back(&lattice_fluid1);
  std::vector<T> constOmegaValues;
  constOmegaValues.push_back(param.omega2);
  constOmegaValues.push_back(param.omega1);
  plint processorLevel = 1;
  integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D<T, DESCRIPTOR>(param.G, constOmegaValues), Box3D(0, param.nx, 0, param.ny, 0, param.nz), blockLattices, processorLevel);

  cout << "Print all the variables\n";
  print_variables();

  pcout << "Reading the geometry file." << endl;

  MultiScalarField3D<int> geometry(param.nx, param.ny, param.nz);
  plb_ifstream geometryFile(param.geometry_fname);

  if (!geometryFile.is_open()) {
      pcout << "Error: could not open geometry file " << fNameIn << endl;
      return -1;
  }
  geometryFile >> geometry;


  cout << "Start simulation\n";

  util::ValueTracer<T> converge1(1.0, lattice_fluid2.getNx(), param.it_con); // it_con parameters velocity/size/threshold
  util::ValueTracer<T> converge2(1.0, lattice_fluid2.getNx(), param.it_con);

  for (plint runs = 1; runs <= param.num_pc; ++runs) { // Loop simulations with varying saturation

      if (runs > 1) {
          pcout << "Using previous Geometry  " << endl;
      }
      else {
          PorousMediaSetup(lattice_fluid1, lattice_fluid2);
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
            << "Starting simulation with rho 1:  " << param.rho_fluid1[runs] << endl;
      pcout << endl
            << "Starting simulation with rho 2:  " << param.rho_fluid2[runs] << endl;

      plint checkconv = 0;
      plint iT = 0;



      while (checkconv == 0) { // Main loop over time iterations.
          iT = iT + 1;

          lattice_fluid1.collideAndStream();
          lattice_fluid2.collideAndStream();

//          if (iT % gifIter == 0) {
//              writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
//              writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
//          }
//          if (iT % vtkIter == 0) {
//              writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
//              writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
//          }

  if (runs==1){

          if (iT == 1) {
//              writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
//              writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);

          }}

converge1.takeValue(getStoredAverageEnergy(lattice_fluid1), true); //check for convergence
converge2.takeValue(getStoredAverageEnergy(lattice_fluid2), true); //check for convergence

          if (iT % param.it_info == 0) {
              pcout << "Iteration:  " << iT << endl;

              if ((converge1.hasConverged()) && (converge2.hasConverged())) {

                  checkconv = 1;
//                  writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
//                  writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
              }


              if ((convergemeanJ1 < convergedmeanJ1) && (convergemeanJ2 < convergedmeanJ2)) {
  //                writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
  //                writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                  pcout << "Simulation has converged" << endl;
  //                writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
  //                writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);

                  //checkconv = 1;
              }
          }

          if (param.it_max == iT) {
              pcout << "Simulation has reached maximum iteration" << endl;
              checkconv = 1;
  //            writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
  //            writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
  //            writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
  //            writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
          }

          initializeAtEquilibrium(lattice_fluid1, Box3D(1, 2, 0, param.ny, 0, param.nz), param.rho_fluid1[runs], Array<T, 3>(0., 0., 0.));
          initializeAtEquilibrium(lattice_fluid2, Box3D(1, 2, 0, param.ny, 0, param.nz), param.rhoNoFluid, Array<T, 3>(0., 0., 0.));
          initializeAtEquilibrium(lattice_fluid1, Box3D(param.nx - 2, param.nx-1, 0, param.ny, 0, param.nz), param.rhoNoFluid, Array<T, 3>(0., 0., 0.));
          initializeAtEquilibrium(lattice_fluid2, Box3D(param.nx - 2, param.nx-1, 0, param.ny, 0, param.nz), param.rho_fluid2[runs], Array<T, 3>(0., 0., 0.));

          lattice_fluid1.initialize();
          lattice_fluid2.initialize();
      }
  }

}




int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outputDir);


    // The try-catch blocks catch exceptions in case an error occurs,
    // and terminate the program properly with a nice error message.

    // 1. Read command-line parameter: the input file name.
    cout << "read command line parameters \n";
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

    // 3. Execute the main program.
    cout << "Executing main program \n";
    try {
        runProgram();
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
