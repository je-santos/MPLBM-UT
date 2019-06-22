/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Wim Degruyter */

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <cstdlib>

using namespace plb;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

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

void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>& geometry, plint run, plint runnum, std::string GeometryName)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();
	plint run_diff = ((runnum - 1)/2)+1;
	std::string fNameIn_temp1 = fNameIn + GeometryName;
	std::string fNameIn_temp = "0";
	
	pcout  <<"\n"   << " Run"<< run << std::endl;
	
    Box3D sliceBox(3,3, 0,ny-1, 0,nz-1);
	
	// lattice_f1_forK_1_.dat
	
	if (run == 1) {  // original geometry - absolute permeability
//	pcout  << "Run geometry "<< std::endl; 
	fNameIn_temp = fNameIn_temp1 + ".dat";
	pcout  << "Run absolute permeability "<< std::endl;
	}
	
	if (run > run_diff) { // Fluid 1 - Krnw
	const plint runner = run - run_diff;
	fNameIn_temp =  fNameIn +"/lattice_f1_forK_" + std::to_string(runner) + "_.dat";	
	pcout  << "Run Krnw "<< std::endl; 
	}
	
	if (run > 1 && run < (run_diff+1)) { // Fluid 2 - Krw
	
	fNameIn_temp =  fNameIn +"/lattice_f2_forK_"+ std::to_string(run-1) + "_.dat";		
	pcout  << "Run Krw "<< std::endl; 
	}
	pcout  << "Geometry name is  "<< fNameIn_temp << std::endl; 
	
    std::auto_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
    plb_ifstream geometryFile(fNameIn_temp.c_str());
	
    for (plint iX=3; iX<nx-4; ++iX) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open the geometry file " << fNameIn_temp << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryFile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX,iX, 0,ny-1, 0,nz-1));
    }

    {
		VtkImageOutput3D<T> vtkOut(createFileName("PorousMedium", run, 6), 1.0);
		
   //     VtkImageOutput3D<T> vtkOut("porousMedium", run, 1.0);
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

void porousMediaSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
        OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
        MultiScalarField3D<int>& geometry, T deltaP)
{
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();
    const plint nz = lattice.getNz();

    pcout << "Definition of inlet/outlet." << std::endl;
    Box3D inlet (15,15, 1,ny-2, 1,nz-2);
    boundaryCondition->addPressureBoundary0N(inlet, lattice);
    setBoundaryDensity(lattice, inlet, (T) 1.);

    Box3D outlet(nx-4,nx-4, 1,ny-2, 1,nz-2);
    boundaryCondition->addPressureBoundary0P(outlet, lattice);
    setBoundaryDensity(lattice, outlet, (T) 1. - deltaP*DESCRIPTOR<T>::invCs2);

 //   pcout << "Definition of the geometry." << std::endl;
    // Where "geometry" evaluates to 1, use bounce-back.
    defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
    // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
    defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

 //   pcout << "Initialization of rho and u." << std::endl;
    initializeAtEquilibrium( lattice, lattice.getBoundingBox(), PressureGradient(deltaP, nx) );

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
            *computeVelocityNorm(lattice, Box3D(4,4, 0,ny-1, 0,nz-1)),
            imSize, imSize );

    // Write velocity-norm at x=nx/2.
    imageWriter.writeScaledGif(createFileName("ux_half", run, 6),
            *computeVelocityNorm(lattice, Box3D(nx/2,nx/2, 0,ny-1, 0,nz-1)),
            imSize, imSize );
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, plint run)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", run, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
}

void computePermeability(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T nu, T deltaP, Box3D domain, T& perm, T& meanU)
{
    pcout << "Computing the permeability." << std::endl;

    // Compute only the x-direction of the velocity (direction of the flow).
    plint xComponent = 0;
    plint nx = lattice.getNx();
	plint ny = lattice.getNy();
    plint nz = lattice.getNz();
	Box3D domain1(15, nx-4, 0, ny-1, 0, nz-1);

    meanU = computeAverage(*computeVelocityComponent(lattice, domain1, xComponent));

    pcout << "Average velocity     = " << meanU                         << std::endl;
    pcout << "Lattice viscosity nu = " << nu                            << std::endl;
    pcout << "Grad P               = " << deltaP/(T)(nx-10)             << std::endl;
	perm = nu*meanU / (deltaP/(T)(nx-20));
  //  pcout << "Permeability         = " << perm 						<< std::endl;
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
	std::string GeometryName ;
	plint maxT;
	
     std::string xmlFname;
  try {
      global::argv(1).read(xmlFname);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: "
              << (std::string)global::argv(0) << " input-file.xml" << std::endl;
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
	
	
	document["folder"]["out_f"].read(fNameOut);	
	document["folder"]["in_f"].read(fNameIn);
	
	document["simulations"]["press"].read(deltaP);
	document["simulations"]["num"].read(Run);
	document["simulations"]["iter"].read(maxT);
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
	pcout << "Total simulations" << runnum << std::endl;

	
	for (plint run = 1; run <= runnum; ++run) {
		
 //   pcout << "Creation of the lattice." << std::endl;
    MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx,ny,nz, new BGKdynamics<T,DESCRIPTOR>(omega));
    // Switch off periodicity.
    lattice.periodicity().toggleAll(false);		

 //   pcout << "Reading the geometry file." << std::endl;
    MultiScalarField3D<int> geometry(nx,ny,nz);
    readGeometry(fNameIn, fNameOut, geometry, run, runnum, GeometryName);

    // pcout << "nu = " << nu << std::endl;
    // pcout << "deltaP = " << deltaP << std::endl;
    // pcout << "omega = " << omega << std::endl;
    // pcout << "nx = " << lattice.getNx() << std::endl;
    // pcout << "ny = " << lattice.getNy() << std::endl;
    // pcout << "nz = " << lattice.getNz() << std::endl;

    porousMediaSetup(lattice, createLocalBoundaryCondition3D<T,DESCRIPTOR>(), geometry, deltaP);

    // The value-tracer is used to stop the simulation once is has converged.
    // 1st parameter:velocity
    // 2nd parameter:size
    // 3rd parameters:threshold
    // 1st and second parameters ae used for the length of the time average (size/velocity)

	
    util::ValueTracer<T> converge(1.0,1000.0,1.0e-5);

    pcout << "Simulation begins" << std::endl;
    plint iT=0;

 //   const plint maxT = 1000;
	
    for (;iT<maxT; ++iT) {
        if (iT % 200 == 0) {
            pcout << "Iteration " << iT << std::endl;
        }
        if (iT % 250 == 0 && iT > 0) {
            writeGifs(lattice,iT,run);
        }

        lattice.collideAndStream();
        converge.takeValue(getStoredAverageEnergy(lattice),true);

        if (converge.hasConverged()) {
            break;
        }
    }

    pcout << "End of simulation at iteration " << iT << " for Run "<< run << std::endl;

 //   pcout << "Permeability:" << std::endl;
    computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(), Perm, Vel);
	
   	
	perm[run]=Perm;
	meanU[run]=Vel;
	
	rel_perm[run]=perm[run]/perm[1];
	pcout << "Relative Permeability = " << rel_perm[run] << std::endl;

    pcout << "Writing VTK file ..." << std::endl;
    writeVTK(lattice, iT, run);
 //   pcout << "Finished!" << std::endl << std::endl;

 //   return 0;
	}
	
	pcout << "Printing outputs" << std::endl;
	std::string outDir = fNameOut + "/";
	std::string output = outDir + GeometryName + ".dat";
	plb_ofstream ofile(output.c_str());
	ofile << "Outputs" << "\n\n";
	ofile << "Krw from run: 2" << "\n" << "Krnw from run: " << (run_diff+1) << std::endl;
	for (plint runs = 1; runs <= runnum; ++runs) {
		
	
	
	ofile << "Run   = " << runs        << std::endl;
	if (runs == 1) {
	ofile << "Absolute Permeability   = " << perm[runs]         << std::endl;
	}
	ofile << "Relative Permeability   = " << rel_perm[runs]         << std::endl;
//	ofile << "Velocity   = " << meanU[runs]  <<"\n"         << std::endl;
	}
}
