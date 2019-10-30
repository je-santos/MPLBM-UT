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

// To do:
//      1) Check if const plint zcomponent = 0 is necessary
//      2) File naming





typedef double T; // Use double-precision arithmetics
// Use a grid which additionally to the f's stores two variables for the external force term.
#define DESCRIPTOR descriptors::ForcedShanChenD3Q19Descriptor

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

  void writeGifs_f2(MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid1,  //creates the pictures
    MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2, plint runs, plint iT)
    {
      const plint imSize = 600;
      const plint zcomponent = 0;
      const plint nx = lattice_fluid2.getNx();
      const plint ny = lattice_fluid2.getNy();
      const plint nz = lattice_fluid2.getNz();
      Box3D slice(0, nx, 0, ny, nz / 2, nz / 2);

      ImageWriter<T> imageWriter("leeloo.map");
      imageWriter.writeScaledGif(createFileName("rho_f2_", 100000000 * runs + iT, 8),
      *computeDensity(lattice_fluid2, slice), imSize, imSize);
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

    T computeVelocity_f1(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1, T nu_f1)
    {
      plint xComponent = 0;
      const plint nx = lattice_fluid1.getNx();
      const plint ny = lattice_fluid1.getNy();
      const plint nz = lattice_fluid1.getNz();

      Box3D domain(3, nx-4, 0, ny-1, 0, nz-1);
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

    void computeRatios(MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid1,
      MultiBlockLattice3D<T,DESCRIPTOR>& lattice_fluid2, T rho_F1, T rho_F2,
      T nu_f1, T nu_f2, T Gads_f1_s2, T Gads_f1_s1, T G, T rhoNoFluid, T meanU1,
      T& meanRho1, T& meanRho2, T& M, T& Ca)
      {
        const plint nx = lattice_fluid1.getNx();
        const plint ny = lattice_fluid1.getNy();
        const plint nz = lattice_fluid1.getNz();

        Box3D domain(3, nx-4, 0, ny-1, 0, nz-1);

        meanRho1 = computeAverageDensity(lattice_fluid1, domain);
        meanRho2 = computeAverageDensity(lattice_fluid2, domain);
        T mu1 = meanRho1*nu_f1;
        T mu2 = meanRho2*nu_f2;

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

      void readGeometry(std::string fNameIn, std::string fNameOut,
        MultiScalarField3D<int>& geometry)
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
          MultiBlockLattice3D<T, DESCRIPTOR>& lattice_fluid2,
          MultiScalarField3D<int>& geometry,
          T rhoNoFluid, T rho_f1, T rho_f2, T Gads_f1_s1, T Gads_f1_s2,
          T Gads_f1_s3, T Gads_f1_s4, T force_f1, T force_f2,
          T nx1_f1, T nx2_f1, T ny1_f1, T ny2_f1, T nz1_f1, T nz2_f1, T nx1_f2,
          T nx2_f2, T ny1_f2, T ny2_f2, T nz1_f2, T nz2_f2, T runs,
          bool load_state, bool print_geom)
          {

            plint nx = lattice_fluid2.getNx();
            plint ny = lattice_fluid2.getNy();
            plint nz = lattice_fluid2.getNz();
            //const T rhoZero = 0.00001;
            //const T deltaP = 0.05;

            pcout << "Definition of the geometry." << endl;
            if (load_state == true) {
              loadBinaryBlock(lattice_fluid1, "lattice_fluid1.dat");
              loadBinaryBlock(lattice_fluid2, "lattice_fluid2.dat");
            }

            if (load_state == false) {
              // First contact angle (labeled with 1)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s1), 1);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s1), 1);

              // NoDynamics (computational efficency, labeled with 2)
              defineDynamics(lattice_fluid1, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);
              defineDynamics(lattice_fluid2, geometry, new NoDynamics<T, DESCRIPTOR>(), 2);

              // Second contact angle (labeled with 3)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s2), 3);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s2), 3);

              // Third contact angle (labeled with 4)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s3), 4);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s3), 4);

              // Fourth contact angle (labeled with 5)
              defineDynamics(lattice_fluid1, geometry, new BounceBack<T, DESCRIPTOR>(-Gads_f1_s4), 5);
              defineDynamics(lattice_fluid2, geometry, new BounceBack<T, DESCRIPTOR>( Gads_f1_s4), 5);
            }


            // Output geometry dynamics
            if (print_geom == true) {
              VtkImageOutput3D<int> vtkOut(createFileName("vtkgeometry", 1, 1), 1.);
              vtkOut.writeData<int>(geometry, "Dynamics", 1.);
              pcout << "Creating geometry vtk file" << endl;
            }

            Array<T, 3> zeroVelocity(0., 0., 0.);

            if (load_state == false) {

              // Initialize  uniform density for target saturation
              pcout << "Initializing Fluids" << endl;

              initializeAtEquilibrium(lattice_fluid2, Box3D(nx1_f2-1, nx2_f2-1,
                                                            ny1_f2-1, ny2_f2-1,
                                                            nz1_f2-1, nz2_f2-1),
                                                            rho_f2, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid1, Box3D(nx1_f2-1, nx2_f2-1,
                                                            ny1_f2-1, ny2_f2-1,
                                                            nz1_f2-1, nz2_f2-1),
                                                            rhoNoFluid, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid1, Box3D(nx1_f1-1, nx2_f1-1,
                                                            ny1_f1-1, ny2_f1-1,
                                                            nz1_f1-1, nz2_f1-1),
                                                            rho_f1, zeroVelocity);

              initializeAtEquilibrium(lattice_fluid2, Box3D(nx1_f1-1, nx2_f1-1,
                                                            ny1_f1-1, ny2_f1-1,
                                                            nz1_f1-1, nz2_f1-1),
                                                            rhoNoFluid, zeroVelocity);
            }

            setExternalVector(lattice_fluid1, lattice_fluid1.getBoundingBox(),
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_f1, 0.));
            setExternalVector(lattice_fluid2, lattice_fluid2.getBoundingBox(),
            DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T, 3>(0., force_f2, 0.));

            lattice_fluid1.initialize();
            lattice_fluid2.initialize();
          }

          int main(int argc, char* argv[])
          {
            // 1. Declaring the variables
            clock_t t;
            t = clock();
            plbInit(&argc, &argv);

            bool load_state ;
            std::string fNameOut ;
            std::string fNameIn  ;
            plint nx ;
            plint ny ;
            plint nz ;

            bool px_f1, py_f1, pz_f1, px_f2, py_f2, pz_f2; //periodicity
            bool pressure_bc;
            plint nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1; //fluid1 configuration
            plint nx1_f2, nx2_f2, ny1_f2, ny2_f2, nz1_f2, nz2_f2; //fluid2 configuration

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

            plint startNum ;

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
            cout << "Reading inputs from xml file \n";
            try {
              XMLreader document(xmlFname);

              document["load_savedstated"].read(load_state);

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

              document["fluids"]["pressure_bc"].read(pressure_bc);

              document["fluids"]["rho_f1_i"].read(rho_f1_inlet);
              document["fluids"]["rho_f2_i"].read(rho_f2_outlet_initial);
              document["fluids"]["rho_f2_f"].read(rho_f2_outlet_final);
              document["fluids"]["rho_d"].read(rhoNoFluid);
              //document["fluids"]["num_pc"].read(rho_f2_step);
              document["fluids"]["drho_f2"].read(drho_f2);


              document["output"]["out_folder"].read(fNameOut);
              document["output"]["save_sim"].read(save_sim);
              document["output"]["convergence"].read(convergence);

              document["output"]["it_max"].read(it_max);
              document["output"]["it_conv"].read(it_conv);
              document["output"]["it_info"].read(it_info);


              document["output"]["it_gif"].read(it_gif);
              document["output"]["it_vtk"].read(it_vtk);
              document["output"]["rho_vtk"].read(rho_vtk);

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
            T M[runnum];
            T Ca[runnum];
            T deltaP[runnum];
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
            T mean_U1[runnum];
            T mean_U2[runnum];
            T mean_rho1[runnum];
            T mean_rho2[runnum];

            const T delta_rho = 0.00005;
            T rho_low = 2 - delta_rho;

            std::string outDir = fNameOut + "/";
            std::string Lattice1 = outDir + fNameIn + "_lattice1.dat";
            std::string Lattice2 = outDir + fNameIn + "_lattice2.dat";

            for (plint readnum = 1; readnum <= runnum; ++readnum) {
              rho_fluid2[readnum] = rho_f2_outlet_initial - (readnum-1)*drho_f2;
              pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
              rho_fluid1[readnum] = rho_f1_inlet;
            }

            const T nu_f1 = ( (T)1 / omega_f1 - 0.5 ) / DESCRIPTOR<T>::invCs2;
            const T nu_f2 = ( (T)1 / omega_f2 - 0.5 ) / DESCRIPTOR<T>::invCs2;

            // Use regularized BGK dynamics to improve numerical stability (but note that BGK dynamics works well too).
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

                  pcout << " The convergence set by the user is = " << convergence << endl;

                  pcout << " The boundary conditions per run are:" << endl;
                  for (plint readnum = 1; readnum <= runnum; ++readnum) {
                    deltaP[readnum]=(rho_fluid1[readnum]-rho_fluid2[readnum])/3;
                    pcout << "Run number = " << readnum << endl;
                    pcout << "Rho_no_1 = " << rho_fluid1[readnum] << endl;
                    pcout << "Rho_no_2 = " << rho_fluid2[readnum] << endl;
                  }

                  pcout << "Reading the geometry file." << endl;
                  MultiScalarField3D<int> geometry(nx, ny, nz);
                  readGeometry(fNameIn, fNameOut, geometry);

                  util::ValueTracer<T> converge1(1.0, lattice_fluid2.getNx(),
                  convergence); // convergence parameters velocity/size/threshold
                  util::ValueTracer<T> converge2(1.0, lattice_fluid2.getNx(),
                  convergence);

                  for (plint runs = 1; runs <= runnum; ++runs) { // Loop simulations with varying saturation
                    pcout << "Run number = " << runs << endl;
                    if (runs > 1)
                    {
                      pcout << "Using previous simulation state  " << endl;
                    }
                    else
                    {
                      PorousMediaSetup(lattice_fluid1, lattice_fluid2, geometry,
                        rhoNoFluid, rho_f1, rho_f2, Gads_f1_s1, Gads_f1_s2,
                        Gads_f1_s3, Gads_f1_s4, force_f1, force_f2,
                        nx1_f1, nx2_f1, ny1_f1, ny2_f1, nz1_f1, nz2_f1,
                        nx1_f2, nx2_f2, ny1_f2, ny2_f2, nz1_f2, nz2_f2, runs,
                        load_state, print_geom);
                      }

                      //T meanJ_2 = 10000.0;
                      //T meanJ_old2;
                      //T meanJ_1 = 10000.0;
                      //T meanJ_old1;


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

                        lattice_fluid1.collideAndStream();
                        lattice_fluid2.collideAndStream();

                        if (iT % it_gif == 0) {
                          writeGifs(  lattice_fluid1, lattice_fluid2, runs + startNum, iT);

                          writeGifs_f2(  lattice_fluid1, lattice_fluid2, runs + startNum, iT);

                          writeGifs2( lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                        }

                        if (iT % it_vtk == 0) {
                          writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
                          writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
                          writeVTK_vel(lattice_fluid1, runs + startNum);
                        }

                        if (iT == 1) {
                          //writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                          //writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                          //writeVTK(lattice_fluid1, runs + startNum,iT, nx, ny, nz);
                          //writeVTK2(lattice_fluid2, runs + startNum,iT, nx, ny, nz);
                        }


                        converge1.takeValue(getStoredAverageDensity(lattice_fluid1), true); //check for convergence
                        converge2.takeValue(getStoredAverageDensity(lattice_fluid2), true); //check for convergence

                        if (iT % it_conv == 0) {
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
                            saveBinaryBlock(lattice_fluid1, Lattice1);
                            saveBinaryBlock(lattice_fluid2, Lattice2);

                            // Calculate velocity here for both fluids in x-direction and pcout

                            T meanU1 = computeVelocity_f1(lattice_fluid1, nu_f1);
                            T meanU2 = computeVelocity_f2(lattice_fluid2, nu_f2);
                            mean_U1[runs] = meanU1;
                            mean_U2[runs] = meanU2;


                            T rho_F1=rho_fluid1[runs];
                            T rho_F2=rho_fluid2[runs];

                            // calculate capillary number & dynamic viscosity ratio
                            computeRatios(lattice_fluid1, lattice_fluid2,
                              rho_F1, rho_F2, nu_f1, nu_f2,
                               Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1,
                               meanRho1, meanRho2, M1, Ca1);

                            M[runs]  = M1;
                            Ca[runs] = Ca1;
                          }

                          if ((convergemeanJ1 < convergedmeanJ1) && (convergemeanJ2 < convergedmeanJ2)) {
                            writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                            writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                            pcout << "Simulation has converged" << endl;
                            writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                            writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);

                            // Calculate velocity here for both fluids in x-direction and pcout
                            T meanU1 = computeVelocity_f1(lattice_fluid1, nu_f1);
                            T meanU2 = computeVelocity_f2(lattice_fluid2, nu_f2);
                            mean_U1[runs] = meanU1;
                            mean_U2[runs] = meanU2;

                            // calculate capillary number & dynamic viscosity ratio
                            computeRatios(lattice_fluid1, lattice_fluid2, rho_F1, rho_F2, nu_f1, nu_f2, Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1, meanRho1, meanRho2, M1, Ca1);
                            M[runs]=M1;
                            Ca[runs]=Ca1;
                          }
                        }

                        if (it_max == iT) {
                          pcout << "Simulation has reached maximum iteration" << endl;
                          checkconv = 1;
                          writeVTK(lattice_fluid1, runs + startNum, iT, nx, ny, nz);
                          writeVTK2(lattice_fluid2, runs + startNum, iT, nx, ny, nz);
                          writeGifs(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                          writeGifs2(lattice_fluid1, lattice_fluid2, runs + startNum, iT);
                          saveBinaryBlock(lattice_fluid1, Lattice1);
                          saveBinaryBlock(lattice_fluid2, Lattice2);

                          // Calculate velocity here for both fluids in x-direction and pcout

                          T meanU1 = computeVelocity_f1(lattice_fluid1, nu_f1);
                          T meanU2 = computeVelocity_f2(lattice_fluid2, nu_f2);
                          mean_U1[runs] = meanU1;
                          mean_U2[runs] = meanU2;


                          T rho_F1=rho_fluid1[runs];
                          T rho_F2=rho_fluid2[runs];

                          // calculate capillary number & dynamic viscosity ratio
                          computeRatios(lattice_fluid1, lattice_fluid2, rho_F1, rho_F2, nu_f1, nu_f2, Gads_f1_s2, Gads_f1_s1, G, rhoNoFluid, meanU1, meanRho1, meanRho2, M1, Ca1);
                          M[runs]=M1;
                          Ca[runs]=Ca1;
                        }

                        // Boundary conditions
                        if (pressure_bc == true)
                        {
                          initializeAtEquilibrium(lattice_fluid1, Box3D(1, 2, 0, ny-1, 0, nz-1), rho_fluid1[runs], Array<T, 3>(0., 0., 0.));
                          initializeAtEquilibrium(lattice_fluid2, Box3D(1, 2, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
                          initializeAtEquilibrium(lattice_fluid1, Box3D(nx - 2, nx-1, 0, ny-1, 0, nz-1), rhoNoFluid, Array<T, 3>(0., 0., 0.));
                          initializeAtEquilibrium(lattice_fluid2, Box3D(nx - 2, nx-1, 0, ny-1, 0, nz-1), rho_fluid2[runs], Array<T, 3>(0., 0., 0.));

                          lattice_fluid1.initialize();
                          lattice_fluid2.initialize();

                        }
                      }
                    }

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
