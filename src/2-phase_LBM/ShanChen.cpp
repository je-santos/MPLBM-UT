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
// Use a grid which additionally to the f's stores two variables for the external force term.
#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor

 //creates the gifs
void writeGif_f1(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
  MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, string runs, plint iT)
  {
    const plint imSize = 600;
    const plint nx = lattice_fluid2.getNx();
    const plint ny = lattice_fluid2.getNy();
    const plint nz = lattice_fluid2.getNz();
    Box3D slice(0, nx, 0, ny, nz / 2, nz / 2);

    string im_name;

    im_name = "rho_f1_";
    im_name.append(runs);
    im_name.append("_");

    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName(im_name, iT, 8),
    *computeDensity(lattice_fluid1, slice), imSize, imSize);
  }


  void writeVTK_vel(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, plint runs)
  {
    VtkImageOutput3D<T> vtkOut(createFileName("vtk_vel_rho1_", 1 * runs, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice_fluid1), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice_fluid1), "velocity", 1.);
  }


  void writeGif_f1_y(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, string runs, plint iT)
    {
      const plint imSize = 600;
      const plint nx = lattice_fluid2.getNx();
      const plint ny = lattice_fluid2.getNy();
      const plint nz = lattice_fluid2.getNz();
      Box3D slice(0, nx, ny / 2, ny / 2, 0, nz);

      string im_name;

      im_name = "rho_f1_y_";
      im_name.append(runs);
      im_name.append("_");


      ImageWriter<T> imageWriter("leeloo.map");
      imageWriter.writeScaledGif(createFileName(im_name, iT, 8),
      *computeDensity(lattice_fluid1, slice), imSize, imSize);
    }

    void writeVTK_rho(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid,
      string im_name, string runs, plint iter, plint nx, plint ny, plint nz)
    {


      im_name.append(runs);
      im_name.append("_");

      const plint zcomponent = 0;
      VtkImageOutput3D<double> vtkOut(createFileName(im_name, iter, 8), 1.);
      vtkOut.writeData<double>((*computeDensity(lattice_fluid)), "Density", 1.);
    }

	 void writeVTK_vel(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid,
      string im_name, string runs, plint iter)
    {
      plint xComponent = 0;
      const plint nx = lattice_fluid.getNx();
      const plint ny = lattice_fluid.getNy();
      const plint nz = lattice_fluid.getNz();
      Box3D domain(0, nx-1, 0, ny-1, 0, nz-1);

      im_name.append(runs);
      im_name.append("_");

      const plint zcomponent = 0;
      VtkImageOutput3D<double> vtkOut(createFileName(im_name, iter, 8), 1.);
      vtkOut.writeData<double>((*computeVelocityComponent(lattice_fluid, domain, xComponent)), "Velocity", 1.);
    }



    T computeVelocity_f1(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, T nu_f1)
    {
      plint xComponent = 0;
      const plint nx = lattice_fluid1.getNx();
      const plint ny = lattice_fluid1.getNy();
      const plint nz = lattice_fluid1.getNz();

      Box3D domain(3, 4, 0, ny-1, 0, nz-1);
      T meanU1 = computeAverage(*computeVelocityComponent(lattice_fluid1, domain, xComponent));
      pcout << "Average velocity for fluid1 in x direction    = "<< meanU1<<std::endl;
      return meanU1;
    }

    T computeVelocity_f2(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid2, T nu_f2)
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

      void readGeometry(std::string fNameIn, std::string fNameOut,
        MultiScalarField3D<int>& geometry)
        {
          const plint nx = geometry.getNx();
          const plint ny = geometry.getNy();
          const plint nz = geometry.getNz();

          Box3D sliceBox(0,0, 0,ny-1, 0,nz-1);
          std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
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
            std::unique_ptr<MultiScalarField3D<T> > floatTags = copyConvert<int,T>(geometry, geometry.getBoundingBox());
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

        void setboundaryvalue(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
                      MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2,
                      Box3D inlet, Box3D outlet,
                      T rho_f1_inlet, T rho_f2_outlet, T rhoNoFluid) {

          setBoundaryDensity(lattice_fluid1, inlet, rho_f1_inlet);
          setBoundaryDensity(lattice_fluid2, inlet, rhoNoFluid);
          setBoundaryDensity(lattice_fluid1, outlet, rhoNoFluid);
          setBoundaryDensity(lattice_fluid2, outlet, rho_f2_outlet);                          

        }

        void PorousMediaSetup(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,
          MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2,
          std::string fluid1_savefile, std::string fluid2_savefile,
          MultiScalarField3D<int>& geometry, 
          OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
          bool bulidin_pressure_bc, Box3D inlet, Box3D outlet,T rho_f1_inlet, T rho_f2_outlet_initial,
          T rhoNoFluid, T rho_f1, T rho_f2, T Gads_f1_s1, T Gads_f1_s2,
          T Gads_f1_s3, T Gads_f1_s4, T force_f1, T force_f2,
          T nx1_f1, T nx2_f1, T ny1_f1, T ny2_f1, T nz1_f1, T nz2_f1, T nx1_f2,
          T nx2_f2, T ny1_f2, T ny2_f2, T nz1_f2, T nz2_f2, T runs,
          bool load_state, bool print_geom)
          {

            plint nx = lattice_fluid2.getNx();
            plint ny = lattice_fluid2.getNy();
            plint nz = lattice_fluid2.getNz();


            pcout << "Definition of the geometry." << endl;

            Array<T, 3> zeroVelocity(0., 0., 0.);
            if (bulidin_pressure_bc==true){
              boundaryCondition->addPressureBoundary0N(inlet, lattice_fluid1);
              boundaryCondition->addPressureBoundary0N(inlet, lattice_fluid2);
              setBoundaryDensity(lattice_fluid1, inlet, rho_f1_inlet);
              setBoundaryDensity(lattice_fluid2, inlet, rhoNoFluid);

              boundaryCondition->addPressureBoundary0P(outlet, lattice_fluid1);
              boundaryCondition->addPressureBoundary0P(outlet, lattice_fluid2);
              setBoundaryDensity(lattice_fluid1, outlet, rhoNoFluid);
              setBoundaryDensity(lattice_fluid2, outlet, rho_f2_outlet_initial);
              delete boundaryCondition;
            }

            if (load_state == true) {

              loadBinaryBlock(lattice_fluid1, fluid1_savefile);
              loadBinaryBlock(lattice_fluid2, fluid2_savefile);

			        plb_ifstream ifile("runnum.dat");
			        if (ifile.is_open()) {
				      ifile >> runs;
			         }
            }

            if (load_state == false) {
              // NoDynamics (computational efficency, labeled with 2)
              defineDynamics(lattice_fluid1, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);
              defineDynamics(lattice_fluid2, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

            }

              // First contact angle (labeled with 1)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s1), 1);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s1), 1);


              // Second contact angle (labeled with 3)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s2), 3);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s2), 3);

              // Mesh contact angle (labeled with 4). Neutral wet
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( 0), 4);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>( 0), 4);

              // Third contact angle (labeled with 4)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s3), 5);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s3), 5);

              // Fourth contact angle (labeled with 5)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s4), 6);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s4), 6);


              if (load_state == false) {
              pcout << "Initializing Fluids" << endl;

              initializeAtEquilibrium(lattice_fluid2, Box3D(nx1_f2-1, nx2_f2-1,
                                                            ny1_f2-1, ny2_f2-1,
                                                            nz1_f2-1, nz2_f2-1),
                                                            rho_f2, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid1, Box3D(nx1_f2-1, nx2_f2-1,
                                                            ny1_f2-1, ny2_f2-1,
                                                            nz1_f2-1, nz2_f2-1),
                                                            rhoNoFluid, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid1, Box3D(nx1_f1, nx2_f1,
                                                            ny1_f1, ny2_f1,
                                                            nz1_f1, nz2_f1),
                                                            rho_f1, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid2, Box3D(nx1_f1, nx2_f1,
                                                            ny1_f1, ny2_f1,
                                                            nz1_f1, nz2_f1),
                                                            rhoNoFluid, zeroVelocity);



            setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(),
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(force_f1, 0., 0.));
            setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(),
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(force_f2, 0., 0.));

            lattice_fluid1.initialize();
            lattice_fluid2.initialize();
          }

            // Output geometry dynamics
            if (print_geom == true) {
              VtkImageOutput3D<int> vtkOut(createFileName("vtkgeometry", 1, 1), 1.);
              vtkOut.writeData<int>(geometry, "Dynamics", 1.);
              pcout << "Creating geometry vtk file" << endl;
            }

          }

          void savestate(std::string outDir,
                  MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1,
                  MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid2,
                  plint iT,plint runs, std::string save_name) {
                  
                  std::string savetime = std::to_string(iT);
                  std::string runnumber = std::to_string(runs);
                  //path to save state
                  std::string Lattice1 = outDir + save_name + "_run_" + runnumber + "_" + savetime + "_lattice1.dat";
                  std::string Lattice2 = outDir + save_name + "_run_" + runnumber + "_" + savetime + "_lattice2.dat";
                  // saves a binary file (heavy) with the sim state
                  saveBinaryBlock(lattice_fluid1, Lattice1);
                  saveBinaryBlock(lattice_fluid2, Lattice2);
          }

          int main(int argc, char* argv[])
          {
            // 1. Declaring the variables
            clock_t t;
            t = clock();
            plbInit(&argc, &argv);

            bool load_state ;
            std::string fluid1_savefile;
            std::string fluid2_savefile;
            std::string fNameOut ;
            std::string fNameIn  ;
            std::string save_name;
            plint nx, ny, nz;

            bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; //periodicity
            bool bulidin_pressure_bc,use_manual_boundary;
            plint nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1; //fluid1 configuration
            plint nx1_f2, nx2_f2, ny1_f2, ny2_f2, nz1_f2, nz2_f2; //fluid2 configuration
            plint inlet_x1, inlet_x2, inlet_y1, inlet_y2, inlet_z1, inlet_z2; 
            plint outlet_x1, outlet_x2, outlet_y1, outlet_y2, outlet_z1, outlet_z2; 
            plint update_frequencey;

            T G ;
            T omega_f1 ;
            T omega_f2 ;
            T force_f1 ;
            T force_f2 ;
            T Gads_f1_s1 ;
            T Gads_f1_s2 ;
            T Gads_f1_s3 ;
            T Gads_f1_s4 ;

            T rho_f1  ;
            T rho_f2  ;

            T rho_f1_inlet ;
            T rho_f2_outlet_initial, rho_f2_outlet_final;
            T rhoNoFluid;
            T rho_f2_step;
            T drho_f2;

            plint it_max ;
            plint it_conv ;
            plint it_info ;
            plint it_vtk ;
            plint it_gif ;
            plint it_save;
            plint restart_run_number;
            plint restart_time;

            bool save_sim, rho_vtk, print_geom, print_stl ;

            T convergence ;

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
            pcout << "Reading inputs from xml file \n";
            try {
              XMLreader document(xmlFname);

              document["restart"]["load_savedstated"].read(load_state);
              document["restart"]["fluid1_savefile"].read(fluid1_savefile);
              document["restart"]["fluid2_savefile"].read(fluid2_savefile);
              document["restart"]["restart_run_number"].read(restart_run_number);
              document["restart"]["restart_time"].read(restart_time);

              document["geometry"]["file_geom"].read(fNameIn);
              document["geometry"]["size"]["x"].read(nx);
              document["geometry"]["size"]["y"].read(ny);
              document["geometry"]["size"]["z"].read(nz);
              document["geometry"]["per"]["fluid1"]["x"].read(px_f1);
              document["geometry"]["per"]["fluid1"]["y"].read(py_f1);
              document["geometry"]["per"]["fluid1"]["z"].read(pz_f1);
              document["geometry"]["per"]["fluid2"]["x"].read(px_f2);
              document["geometry"]["per"]["fluid2"]["y"].read(py_f2);
              document["geometry"]["per"]["fluid2"]["z"].read(pz_f2);

              document["init"]["fluid1"]["x1"].read(nx1_f1);
              document["init"]["fluid1"]["x2"].read(nx2_f1);
              document["init"]["fluid1"]["y1"].read(ny1_f1);
              document["init"]["fluid1"]["y2"].read(ny2_f1);
              document["init"]["fluid1"]["z1"].read(nz1_f1);
              document["init"]["fluid1"]["z2"].read(nz2_f1);

              document["init"]["fluid2"]["x1"].read(nx1_f2);
              document["init"]["fluid2"]["x2"].read(nx2_f2);
              document["init"]["fluid2"]["y1"].read(ny1_f2);
              document["init"]["fluid2"]["y2"].read(ny2_f2);
              document["init"]["fluid2"]["z1"].read(nz1_f2);
              document["init"]["fluid2"]["z2"].read(nz2_f2);

              document["fluids"]["Gc"].read(G);
              document["fluids"]["omega_f1"].read(omega_f1);
              document["fluids"]["omega_f2"].read(omega_f2);

              document["fluids"]["force_f1"].read(force_f1);
              document["fluids"]["force_f2"].read(force_f2);

              document["fluids"]["G_ads_f1_s1"].read(Gads_f1_s1);
              document["fluids"]["G_ads_f1_s2"].read(Gads_f1_s2);
              document["fluids"]["G_ads_f1_s3"].read(Gads_f1_s3);
              document["fluids"]["G_ads_f1_s4"].read(Gads_f1_s4);

              document["fluids"]["rho_f1"].read(rho_f1);
              document["fluids"]["rho_f2"].read(rho_f2);

              document["fluids"]["bulidin_pressure_bc"].read(bulidin_pressure_bc);
              document["fluids"]["manual_boundary"]["use_manual_boundary"].read(use_manual_boundary);
              document["fluids"]["manual_boundary"]["update_frequence"].read(update_frequencey);
              document["fluids"]["inlet_position"]["x1"].read(inlet_x1);
              document["fluids"]["inlet_position"]["x2"].read(inlet_x2);
              document["fluids"]["inlet_position"]["y1"].read(inlet_y1);
              document["fluids"]["inlet_position"]["y2"].read(inlet_y2);
              document["fluids"]["inlet_position"]["z1"].read(inlet_z1);
              document["fluids"]["inlet_position"]["z2"].read(inlet_z2);
              document["fluids"]["outlet_position"]["x1"].read(outlet_x1);
              document["fluids"]["outlet_position"]["x2"].read(outlet_x2);
              document["fluids"]["outlet_position"]["y1"].read(outlet_y1);
              document["fluids"]["outlet_position"]["y2"].read(outlet_y2);
              document["fluids"]["outlet_position"]["z1"].read(outlet_z1);
              document["fluids"]["outlet_position"]["z2"].read(outlet_z2);

              document["fluids"]["rho_f1_i"].read(rho_f1_inlet);
              document["fluids"]["rho_f2_i"].read(rho_f2_outlet_initial);
              document["fluids"]["rho_f2_f"].read(rho_f2_outlet_final);
              document["fluids"]["rho_d"].read(rhoNoFluid);
              document["fluids"]["drho_f2"].read(drho_f2);


              document["output"]["out_folder"].read(fNameOut);
              document["output"]["convergence"].read(convergence);

              document["output"]["it_max"].read(it_max);
              document["output"]["it_conv"].read(it_conv);

              document["output"]["it_gif"].read(it_gif);
              document["output"]["it_vtk"].read(it_vtk);
              document["output"]["rho_vtk"].read(rho_vtk);
              document["output"]["it_save"].read(it_save);
              document["output"]["save_sim"].read(save_sim);
              document["output"]["save_name"].read(save_name);

              document["output"]["print_geom"].read(print_geom);
              document["output"]["print_stl"].read(print_stl);

            }
            catch (PlbIOException& exception) {
              pcout << exception.what() << std::endl;
              pcout << exception.what() << std::endl;
              return -1;
            }

            plint runnum = ((rho_f2_outlet_initial - rho_f2_outlet_final)/drho_f2)+1;
            global::directories().setOutputDir(fNameOut);

            T rho_fluid1[runnum] ;
            T rho_fluid2[runnum] ;
            T deltaP[runnum];
            T k1_high;
            T k2_high;
            T meanRho1;
            T meanRho2;
            T mu1;
            T mu2;
            T rho_F1;
            T rho_F2;
            T mean_U1[runnum];
            T mean_U2[runnum];
            T mean_rho1[runnum];
            T mean_rho2[runnum];

            std::string outDir   =   fNameOut;
            std::string Lattice1 =   fNameOut + "lattice1.dat";
            std::string Lattice2 =   fNameOut + "lattice2.dat";

            for (plint readnum = 1; readnum <= runnum; ++readnum) {
              rho_fluid2[readnum] = rho_f2_outlet_initial - (readnum-1)*drho_f2;
              pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
              rho_fluid1[readnum] = rho_f1_inlet;
            }

            rho_fluid2[1]=rho_fluid2[1]+0.02; // Adding 0.02 to fluid 2 at outlet in first run to prevent premature fluid invasion

            const T nu_f1 = ( (T)1 / omega_f1 - 0.5 ) / DESCRIPTOR<T>::invCs2;
            const T nu_f2 = ( (T)1 / omega_f2 - 0.5 ) / DESCRIPTOR<T>::invCs2;

            // Use regularized BGK dynamics to improve numerical stability
            // (but note that BGK dynamics works well too).
            MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid2( nx, ny, nz,
              new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega_f2));
            MultiBlockLattice3D<T, DESCRIPTOR> lattice_fluid1( nx, ny, nz,
              new ExternalMomentRegularizedBGKdynamics<T, DESCRIPTOR>(omega_f1));

                lattice_fluid2.periodicity().toggle(0, px_f2);
                lattice_fluid1.periodicity().toggle(0, px_f1);
                lattice_fluid2.periodicity().toggle(1, py_f2);
                lattice_fluid1.periodicity().toggle(1, py_f1);
                lattice_fluid2.periodicity().toggle(2, pz_f2);
                lattice_fluid1.periodicity().toggle(2, pz_f1);

                vector<MultiBlockLattice3D<T, DESCRIPTOR>*> blockLattices;
                blockLattices.push_back( &lattice_fluid2 );
                blockLattices.push_back( &lattice_fluid1 );

                std::vector<T> constOmegaValues;
                constOmegaValues.push_back( omega_f2 );
                constOmegaValues.push_back( omega_f1 );
                plint processorLevel = 1;

                integrateProcessingFunctional(new ShanChenMultiComponentProcessor3D<T,
                DESCRIPTOR>(G, constOmegaValues), Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                blockLattices, processorLevel);

                  pcout << "The convergence set by the user is = " << convergence << endl;

                  if (bulidin_pressure_bc == true || use_manual_boundary == true) {
                    pcout << "The boundary conditions per run are:" << endl;
                    for (plint readnum = 1; readnum <= runnum; ++readnum) {
                      deltaP[readnum]=(rho_fluid1[readnum]-rho_fluid2[readnum])/3;
                      pcout << "Run number = " << readnum << endl;
                      pcout << "Rho_no_1 = " << rho_fluid1[readnum] << endl;
                      pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
                    }
                }

                  pcout << "Reading the geometry file." << endl;
                  MultiScalarField3D<int> geometry(nx, ny, nz);
                  readGeometry(fNameIn, fNameOut, geometry);

                  //  iteration time steps
                  plint iT = 0;
                  // logging the current run type to see if it's the first run of a restart simulation
                  bool restart_sim = false;
                  //TODO RESTORE BOUNDARY VALUE WHEN LOADING A STATE
                  if (load_state == true) {
                    iT = restart_time;
                    restart_sim == true;
                    pcout << "simulation restart at iT = " << iT <<"\n"<< std::endl;
                    pcout <<"simulation restart at run number: " << restart_run_number <<"\n"<< std::endl;
                  }
                  else 
                  {
                    //set restart_run_number to 1 if it's not a restart simulation
                    restart_run_number = 1;
                  }

                  Box3D inlet(inlet_x1-1, inlet_x2-1, inlet_y1-1, inlet_y2-1, inlet_z1-1, inlet_z2-1);
                  Box3D outlet(outlet_x1-1, outlet_x2-1, outlet_y1-1, outlet_y2-1, outlet_z1-1, outlet_z2-1);
                  // Loop simulations with varying saturation
                  for (plint runs = restart_run_number; runs <= runnum; ++runs) {

                    // turn off stats for efficency
                    lattice_fluid1.toggleInternalStatistics(false);
                    lattice_fluid2.toggleInternalStatistics(false);

                    // save a str for figure naming
                    stringstream save_str;
                    save_str << std::setw(3) << std::setfill('0') << runs;
                    string runs_str;
                    save_str  >> runs_str;

                    pcout << "Run number = " << runs << endl;

                    // re-use the final state of the previous run
                    // if it is a restart run it should setup the domain first
                    if (runs > 1 && restart_sim == false)
                    {
                      pcout << "Using previous simulation state  " << endl;
                      if (!bulidin_pressure_bc){
                          pcout << "Update boundary value " << endl;
                          setboundaryvalue(lattice_fluid1,lattice_fluid2,inlet,outlet,rho_fluid1[runs],rho_fluid2[runs],rhoNoFluid);
                      }
                    }
                    // set-up a new simulation domain
                    else
                    {
                      PorousMediaSetup(lattice_fluid1, lattice_fluid2, fluid1_savefile,
                        fluid2_savefile,geometry,  createLocalBoundaryCondition3D<T,DESCRIPTOR>(), bulidin_pressure_bc, inlet, outlet,
                        rho_f1_inlet, rho_f2_outlet_initial, rhoNoFluid, rho_f1, rho_f2, Gads_f1_s1, Gads_f1_s2,
                        Gads_f1_s3, Gads_f1_s4, force_f1, force_f2,
                        nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1,
                        nx1_f2, nx2_f2, ny1_f2, ny2_f2, nz1_f2, nz2_f2, runs,
                        load_state, print_geom);
                        // PorousMediaSetup is done for a restart sim
                        restart_sim = false;
                      }

                      T new_avg_f1, new_avg_f2, old_avg_f1, old_avg_f2;
                      T relE_f1, relE_f2;

                      pcout << endl
                      << "Starting simulation with rho 1:  " << rho_fluid1[runs] << endl;
                      pcout << endl
                      << "Starting simulation with rho 2:  " << rho_fluid2[runs] << endl;

                      plint checkconv = 0;

                      while (checkconv == 0) { // Main loop over time iterations.
                        iT = iT + 1;

                        // turn on stats to check convergence
                        if (iT % it_conv == 0) {
                          lattice_fluid1.toggleInternalStatistics(true);
                          lattice_fluid2.toggleInternalStatistics(true);
                          }

                        lattice_fluid1.collideAndStream();
                        lattice_fluid2.collideAndStream();

                        // save gifs
                        if (iT % it_gif == 0) {
                          writeGif_f1(  lattice_fluid1, lattice_fluid2, runs_str, iT);
                          writeGif_f1_y(lattice_fluid1, lattice_fluid2, runs_str, iT);
                        }

                        // save vtks
                        if (iT % it_vtk == 0) {
                          writeVTK_rho(lattice_fluid1, "rho_f1_", runs_str, iT, nx, ny, nz);
                          if (rho_vtk == true)
                          {
                            writeVTK_rho(lattice_fluid2, "rho_f2_", runs_str, iT, nx, ny, nz);
                          }
                        }
                        if (save_sim == true) {
                          if (iT % it_save == 0 || iT == it_max) {
                            // saves a binary file (heavy) with the sim state
                            savestate(outDir,lattice_fluid1, lattice_fluid2, iT, runs, save_name);
                          }
                        }
                        if (iT % it_conv == 0 ) {
                          // calculate average change in mass if bcs == pressure
                          new_avg_f1 = getStoredAverageDensity(lattice_fluid1)*(nx*ny*nz);
                          new_avg_f2 = getStoredAverageDensity(lattice_fluid2)*(nx*ny*nz);

                          if (bulidin_pressure_bc == false && use_manual_boundary == false){
                          // calculate average change in momentum if bcs == force
                            new_avg_f1 = getStoredAverageEnergy(lattice_fluid1);
                            new_avg_f2 = getStoredAverageEnergy(lattice_fluid2);
                          }

                          //mean_rho1[runs] = getStoredAverageDensity<T>(lattice_fluid1);
                          //mean_rho2[runs] = getStoredAverageDensity<T>(lattice_fluid2);

                          lattice_fluid1.toggleInternalStatistics(false);
                          lattice_fluid2.toggleInternalStatistics(false);

                          // calculate relative difference
                          relE_f1 = std::fabs(old_avg_f1-new_avg_f1)*100/old_avg_f1;
                          relE_f2 = std::fabs(old_avg_f2-new_avg_f2)*100/old_avg_f2;

                          pcout << "Run num " << runs;
                          pcout << ", Iteration " << iT << std::endl;
                          pcout << "-----------------"  << std::endl;
                          pcout << "Relative difference fluid1: " << setprecision(3)
                          << relE_f1 <<" %"<<std::endl;
                          pcout << "Relative difference fluid2: " << setprecision(3)
                          << relE_f2 <<" %"<<std::endl;
                          pcout << "Has fluid 1 converged?: "<< ((relE_f1 < convergence) ? "TRUE" : "FALSE") <<  std::endl;
                          pcout << "Has fluid 2 converged?: "<< ((relE_f2 < convergence) ? "TRUE" : "FALSE") <<  std::endl;
                          pcout << "-----------------"  << std::endl;

                          // store new properties
                          old_avg_f1 = new_avg_f1;
                          old_avg_f2 = new_avg_f2;

                          if ( relE_f1 < convergence && relE_f2 < convergence ){
                            checkconv = 1;
                            pcout << "Pressure increment has converged" << endl;
                          }
                        }

                        if (it_max == iT) {
                          pcout << "Simulation has reached maximum iteration" << endl;
                          checkconv = 1;
                        }

                          if ( checkconv == 1 ) {
                            writeGif_f1(  lattice_fluid1, lattice_fluid2, runs_str, iT);
                            writeGif_f1_y(lattice_fluid1, lattice_fluid2, runs_str, iT);

                          // saves converged state vtks
                          if ( it_vtk < 100000 ) {

                            writeVTK_rho(lattice_fluid1, "rho_f1_", runs_str, iT, nx, ny, nz);
							              //writeVTK_vel(lattice_fluid1, "vel_f1_", runs_str, iT);

                            if (rho_vtk == true)
                            {
                              writeVTK_rho(lattice_fluid2, "rho_f2_", runs_str, iT, nx, ny, nz);
                            }

                          }

                            // saves a .dat file with the run number (for restarting sim)
                            string run_name;
                            run_name = outDir + "/runnum.dat";
                            plb_ofstream ofile1( run_name.c_str()  );
              							ofile1 << runs << endl;

                            // saves a .dat file (lightweight) with the density
                            string rho_name;
                            rho_name = outDir + "/rho_f1_" + runs_str + ".dat";
  					                plb_ofstream ofile2( rho_name.c_str() );
                            ofile2 << setprecision(2) <<*computeDensity(lattice_fluid1) << endl;

                            // saves a .dat file (lightweight) with the velocity
                            string vel_name;
                            vel_name = outDir + "/vel_f1_" + runs_str + ".dat";
                            plb_ofstream ofile3( vel_name.c_str() );
                            ofile3 << setprecision(1) <<*computeVelocity(lattice_fluid1) << endl;
                            
                            // saves a binary file (heavy) with the sim state
                            if (save_sim == true) {      

                              string converged_save_name = save_name + "_converged";
                              savestate(outDir,lattice_fluid1, lattice_fluid2, iT, runs, converged_save_name);

                            }

                            // Calculate velocity here for both fluids in x-direction
                            T meanU1 = computeVelocity_f1(lattice_fluid1, nu_f1);
                            T meanU2 = computeVelocity_f2(lattice_fluid2, nu_f2);
                            mean_U1[runs] = meanU1;
                            mean_U2[runs] = meanU2;

                            T rho_F1=rho_fluid1[runs];
                            T rho_F2=rho_fluid2[runs];

                          }


                        if (use_manual_boundary == true && iT % update_frequencey == 0) // Pressure boundary conditions
                        {

                          Array<T, 3> zeroVelocity(0., 0., 0.);
                          initializeAtEquilibrium(lattice_fluid1, inlet, rho_fluid1[runs], zeroVelocity);
                          initializeAtEquilibrium(lattice_fluid2, inlet, rhoNoFluid, zeroVelocity);
                          initializeAtEquilibrium(lattice_fluid1, outlet, rhoNoFluid, zeroVelocity);
                          initializeAtEquilibrium(lattice_fluid2, outlet, rho_fluid2[runs], zeroVelocity);
                          
                          lattice_fluid1.initialize();
                          lattice_fluid2.initialize();

                        }
                      }
                    }

                    std::string output = outDir  + "/output.dat";
                    t = clock() - t;
                    pcout << "Simulation took seconds:" << ((float)t)/CLOCKS_PER_SEC << std::endl;
                    plb_ofstream ofile(output.c_str());
                    ofile << "Output of the Simulation Run" << "\n\n";
                    ofile << "Simulation took seconds =" << ((float)t)/CLOCKS_PER_SEC <<"\n" << endl;

          					ofile << "Kinematic viscosity f1 = " << nu_f1 <<"\n" << endl;
          					ofile << "Kinematic viscosity f2 = " << nu_f2 <<"\n" << endl;
          					ofile << "Gads_f1_s1 = " << Gads_f1_s1 <<"\n" << endl;
          					ofile << "Gads_f1_s2 = " << Gads_f1_s2 <<"\n" << endl;
          					ofile << "Gc = " << G <<"\n" << endl;
          					ofile << "Dissolved density = " << rhoNoFluid <<"\n" << endl;
          					ofile << "Inlet density = " << rho_f1_inlet <<"\n" << endl;
          					ofile << "Geometry flow length = " << nx <<"\n" << endl;

                    for (plint runs = 1; runs <= runnum; ++runs) {

                      pcout << "Run    = " << runs       << std::endl;
                      pcout << "Pressure difference =  " << deltaP[runs] << std::endl;

                      ofile << "Run = " << runs << "\n" << endl;
                      ofile << "Pressure difference = " << deltaP[runs] <<"\n" << endl;

                    }
                  }
