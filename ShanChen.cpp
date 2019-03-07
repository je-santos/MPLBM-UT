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



// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param
{

    T omega1, omega2;
    T nu1,nu2;                               // Kinematic viscosity.
    //T lx, ly, lz;                       // Size of computational domain, in physical units.



    T it_max, it_con,it_info, it_gif, it_vtk;
    //plint resolution;                   // Number of lattice nodes along a reference length.
    plint nx, ny, nz;                   // Grid resolution of bounding box.
    plint nx1f1,nx2f1,ny1f1,ny2f1,nz1f1,nz2f1;
    plint nx1f2,nx2f2,ny1f2,ny2f2,nz1f2,nz2f2;

    plint num_pc;
    T rho1_i,rho2_i, rho1_f,rho2_f, rho1,rho2, rho_d;

    T G;
    T force_fluid1,force_fluid2;
    T Gads_f1_s1,Gads_f1_s2;
    bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2;
    std::string input_file, output_folder;



    //T omega;                            // Relaxation parameter.
    //T dx, dt;                           // Discrete space and time steps.
    //plint maxIter, statIter;            // Time for events in lattice units.
    //plint imageIter, vtkIter;
    //bool useParticles;                  // Simulate particles or not.
    //int particleTimeFactor;             // If the particle time factor is 2, then the integration time step
                                        //   for the particles is twice that of the fluid.
    //T particleProbabilityPerCell;       // Probability of injection of a particle at an injection cell at each time step.
    //T cutOffSpeedSqr;                   // Criterion to eliminate particles with very small velocity.
    //int maxNumParticlesToWrite;         // Maximum number of particles in the output VTK files.

    //T outletSpongeZoneWidth;            // Width of the outlet sponge zone.
    //plint numOutletSpongeCells;         // Number of the lattice nodes contained in the outlet sponge zone.
    //int outletSpongeZoneType;           // Type of the outlet sponge zone (Viscosity or Smagorinsky).
    //T targetSpongeCSmago;               // Target Smagorinsky parameter at the end of the Smagorinsky sponge Zone.
    //plint initialIter;                  // Number of initial iterations until the inlet velocity reaches its final value.

    //Box3D inlet, outlet, lateral1;      // Outer domain boundaries in lattice units.
    //Box3D lateral2, lateral3, lateral4;

    //std::string file_geom;
    //std::string folder_output;

    Param()
    { }

    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);

        document["geometry"]["file_geom"].read(input_file);
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

        //px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2;

/*


        document["geometry"]["filename"].read(geometry_fname);
        document["geometry"]["center"]["x"].read(cx);
        document["geometry"]["center"]["y"].read(cy);
        document["geometry"]["center"]["z"].read(cz);
        document["geometry"]["freeSlipWall"].read(freeSlipWall);
        document["geometry"]["lateralFreeSlip"].read(lateralFreeSlip);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

        document["numerics"]["nu"].read(nu);
        document["numerics"]["inletVelocity"].read(inletVelocity);
        document["numerics"]["resolution"].read(resolution);
        document["numerics"]["uLB"].read(uLB);
        document["numerics"]["useSmago"].read(useSmago);
        if (useSmago) {
            document["numerics"]["cSmago"].read(cSmago);
        }

        document["numerics"]["useParticles"].read(useParticles);
        if (useParticles) {
            document["numerics"]["particleTimeFactor"].read(particleTimeFactor);
            document["numerics"]["particleProbabilityPerCell"].read(particleProbabilityPerCell);
            document["numerics"]["cutOffSpeedSqr"].read(cutOffSpeedSqr);
            document["numerics"]["maxNumParticlesToWrite"].read(maxNumParticlesToWrite);
        }

        document["numerics"]["outletSpongeZoneWidth"].read(outletSpongeZoneWidth);
        std::string zoneType;
        document["numerics"]["outletSpongeZoneType"].read(zoneType);
        if ((util::tolower(zoneType)).compare("viscosity") == 0) {
            outletSpongeZoneType = 0;
        } else if ((util::tolower(zoneType)).compare("smagorinsky") == 0) {
            outletSpongeZoneType = 1;
        } else {
            pcout << "The sponge zone type must be either \"Viscosity\" or \"Smagorinsky\"." << std::endl;
            exit(-1);
        }
        document["numerics"]["targetSpongeCSmago"].read(targetSpongeCSmago);

        document["numerics"]["initialIter"].read(initialIter);

        document["output"]["maxT"].read(maxT);
        document["output"]["statT"].read(statT);
        document["output"]["imageT"].read(imageT);
        document["output"]["vtkT"].read(vtkT);
*/

        computeLBparameters();
    }

    void computeLBparameters()
    {

    T rho_fluid1[num_pc];
    T rho_fluid2[num_pc];
    rho_fluid1[0]=rho1_i;

    rho_fluid1[num_pc-1]=rho1_f;
    for (plint i=1; i<num_pc; i++){
      rho_fluid1[i]=rho1_i+(rho1_f-rho1_i)/(num_pc-1)*i;
    }
    const T nu1 = ((T)1 / omega1 - 0.5) / DESCRIPTOR<T>::invCs2;
    const T nu2 = ((T)1 / omega2 - 0.5) / DESCRIPTOR<T>::invCs2;

    //Should this be on a separate init function? We'll see

    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega2));
    MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1(nx, ny, nz, new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega1));

    lattice_fluid1.periodicity().toggle(0, px_f1);
    lattice_fluid2.periodicity().toggle(0, px_f2);
    lattice_fluid1.periodicity().toggle(1, py_f1);
    lattice_fluid2.periodicity().toggle(1, py_f2);
    lattice_fluid1.periodicity().toggle(2, pz_f1);
    lattice_fluid2.periodicity().toggle(2, pz_f2);

    vector<MultiBlockLattice3D<T, DESCRIPTOR>*> blockLattices;
    blockLattices.push_back(&lattice_fluid2);
    blockLattices.push_back(&lattice_fluid1);
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega2);
    constOmegaValues.push_back(omega1);
    plint processorLevel = 1;
    integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D<T, DESCRIPTOR>(G, constOmegaValues), Box3D(0, nx, 0, ny, 0, nz), blockLattices, processorLevel);
    MultiScalarField3D<int> geometry(nx, ny, nz);
    plb_ifstream geometryFile(input_file.c_str());

    if (!geometryFile.is_open()) {
        pcout << "Fatal error: could not open geometry file " << input_file << endl;
        //return -1;
    }
    geometryFile >> geometry;

    for (plint it_dens=0; it_dens<num_pc; ++it_dens){ //main loop
      cout<< it_dens<<"Here we go..." << endl;
      if (it_dens > 0) {
          pcout << "Using fluid distribution  " << endl;
      }
      else{

         PorousMediaSetup(lattice_fluid1,lattice_fluid2,geometry,T rho_d,
          T rho_fluid1, T rho_fluid2, T Gads_f1_s1, T Gads_f1_s2, T force_fluid1,
           T force_fluid2, T nx1f1, T nx2f1, T ny1f1, T ny2f1, T nz1f1, T nz2f1, T nx1f2,
            T nx2f2, T ny1f2, T ny2f2, T nz1f2, T nz2f2));
      }

      cout<< it_dens<< endl;
    }


/*

        dx = ly / (resolution - 1.0);
        dt = (uLB/inletVelocity) * dx;
        T nuLB = nu * dt/(dx*dx);
        omega = 1.0/(DESCRIPTOR<T>::invCs2*nuLB+0.5);
        if (lateralFreeSlip) {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx) + 1;
            nz = util::roundToInt(lz/dx) + 1;
        } else {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx);
            nz = util::roundToInt(lz/dx);
        }
        cxLB = util::roundToInt(cx/dx);
        cyLB = util::roundToInt(cy/dx);
        czLB = util::roundToInt(cz/dx);
        maxIter   = util::roundToInt(maxT/dt);
        statIter  = util::roundToInt(statT/dt);
        imageIter = util::roundToInt(imageT/dt);
        vtkIter   = util::roundToInt(vtkT/dt);
        numOutletSpongeCells = util::roundToInt(outletSpongeZoneWidth/dx);

        inlet    = Box3D(0,      0,      0,      ny-1,   0,      nz-1);
        outlet   = Box3D(nx-1,   nx-1,   0,      ny-1,   0,      nz-1);
        lateral1 = Box3D(1,      nx-2,   0,      0,      0,      nz-1);
        lateral2 = Box3D(1,      nx-2,   ny-1,   ny-1,   0,      nz-1);
        lateral3 = Box3D(1,      nx-2,   1,      ny-2,   0,      0);
        lateral4 = Box3D(1,      nx-2,   1,      ny-2,   nz-1,   nz-1);
    }
*/
}
};

Param param;



// Write VTK file for the flow around the obstacle, to be viewed with Paraview.
//void writeVTK(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& bc, plint iT)
//{
//    VtkImageOutput3D<T> vtkOut(createFileName("volume", iT, PADDING));
//    vtkOut.writeData<float>( *bc.computeVelocityNorm(param.boundingBox()),
//                             "velocityNorm", param.dx/param.dt );
//    vtkOut.writeData<3,float>(*bc.computeVelocity(param.boundingBox()), "velocity", param.dx/param.dt);
//    vtkOut.writeData<float>( *bc.computePressure(param.boundingBox()),
//                             "pressure", param.dx*param.dx/(param.dt*param.dt) );
//}

// Write PPM images on slices.
//void writePPM(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& bc, plint iT)
//{
//    Box3D xSlice(param.cxLB, param.cxLB, 0,          param.ny-1, 0,          param.nz-1);
//    Box3D ySlice(0,          param.nx-1, param.cyLB, param.cyLB, 0,          param.nz-1);
//    Box3D zSlice(0,          param.nx-1, 0,          param.ny-1, param.czLB, param.czLB);

//    ImageWriter<T> writer("leeloo");
//    writer.writeScaledPpm(createFileName("vnorm_xslice", iT, PADDING), *bc.computeVelocityNorm(xSlice));
//    writer.writeScaledPpm(createFileName("vnorm_yslice", iT, PADDING), *bc.computeVelocityNorm(ySlice));
//    writer.writeScaledPpm(createFileName("vnorm_zslice", iT, PADDING), *bc.computeVelocityNorm(zSlice));
//}

void PorousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, MultiScalarField3D<int>& geometry, struct param)
    {
      plint nx = lattice_fluid2.getNx();
      plint ny = lattice_fluid2.getNy();
      plint nz = lattice_fluid2.getNz();

      pcout << "Definition of the geometry." << endl;
      for (plint itX = 0; itX < nx; ++itX) {
          for (plint itY = 0; itY < ny; ++itY) {
              for (plint itZ = 0; itZ < nz; ++itZ) {
                  if (geometry.get(itX, itY, itZ) == 1) {
                      defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(Gads_f1_s1));
                      defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new BounceBack<T, DESCRIPTOR>(-Gads_f1_s1));
                  }
                  if (geometry.get(itX, itY, itZ) == 2) {
                      defineDynamics(lattice_fluid1, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
                      defineDynamics(lattice_fluid2, Box3D(itX, itX, itY, itY, itZ, itZ), new NoDynamics<T, DESCRIPTOR>());
                  }
                  if (geometry.get(itX, itY, itZ) == 3) {
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



void runProgram()
{
  cout<<"Prog here";

}



int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(param.output_folder);


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
    try {
        param = Param(xmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }

    // 3. Execute the main program.
    try {
        runProgram();
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}
