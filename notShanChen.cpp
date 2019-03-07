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

static std::string outputDir("./tmp/");

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param
{
    T nu;                               // Kinematic viscosity.
    T lx, ly, lz;                       // Size of computational domain, in physical units.
    T cx, cy, cz;                       // Position of the center of the obstacle, in physical units.
    plint cxLB, cyLB, czLB;             // Position of the center of the obstacle, in lattice units.
    bool freeSlipWall;                  // Use free-slip condition on obstacle, as opposed to no-slip?
    bool lateralFreeSlip;               // Use free-slip lateral boundaries or periodic ones?
    T maxT, statT, imageT, vtkT;        // Time, in physical units, at which events occur.
    plint resolution;                   // Number of lattice nodes along a reference length.
    T inletVelocity;                    // Inlet x-velocity in physical units.
    T uLB;                              // Velocity in lattice units (numerical parameters).
    bool useSmago;                      // Use a Smagorinsky LES model or not.
    T cSmago;                           // Parameter for the Smagorinsky LES model.
    plint nx, ny, nz;                   // Grid resolution of bounding box.
    T omega;                            // Relaxation parameter.
    T dx, dt;                           // Discrete space and time steps.
    plint maxIter, statIter;            // Time for events in lattice units.
    plint imageIter, vtkIter;
    bool useParticles;                  // Simulate particles or not.
    int particleTimeFactor;             // If the particle time factor is 2, then the integration time step
                                        //   for the particles is twice that of the fluid.
    T particleProbabilityPerCell;       // Probability of injection of a particle at an injection cell at each time step.
    T cutOffSpeedSqr;                   // Criterion to eliminate particles with very small velocity.
    int maxNumParticlesToWrite;         // Maximum number of particles in the output VTK files.

    T outletSpongeZoneWidth;            // Width of the outlet sponge zone.
    plint numOutletSpongeCells;         // Number of the lattice nodes contained in the outlet sponge zone.
    int outletSpongeZoneType;           // Type of the outlet sponge zone (Viscosity or Smagorinsky).
    T targetSpongeCSmago;               // Target Smagorinsky parameter at the end of the Smagorinsky sponge Zone.
    plint initialIter;                  // Number of initial iterations until the inlet velocity reaches its final value.

    Box3D inlet, outlet, lateral1;      // Outer domain boundaries in lattice units.
    Box3D lateral2, lateral3, lateral4;

    std::string geometry_fname;

    Param()
    { }

    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);
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

        //computeLBparameters();
    }

    void computeLBparameters()
    {
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

void runProgram()
{
  cout<<"Prog here";
  cout<<param.cx;
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
