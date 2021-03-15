// Last update: December 1, 2017
// This is main for the RDF computation of optionally polydisperse systems.

//  2017.12.18 NB Notes:  Program Outline
/*
 *  This code has been shown to give identical results to the single-timestep computation done by Ovito as of Dec 2017.
 *
 *  Input (in real units, Molar) the concentrations and number of each particle type.
 *    It then calculates the box length based on the above quantities as specified for type-1.  Verify it's correct!
 *
 *  The code requires instantaneous dump files (1 timestep per file). [See "Necessary Code Improvements" section below.]
 *    Flags specify if instantaneous files are to be generated from a larger movie file &/or deleted after generation.
 *      Movie file must have its name specified internally in the 'generate_inst_dump_files(~)' command.
 *        [Functionality to read directly from movie file without intermediate step desirable, see "Necessary Code Improvemenets" comment section below.]
 *
 *  It will use (and if flagged, generate) 'dataSetCount' instantaneous dump files, starting with the input 'initDumpStep'.
 *
 *      e.g. if 'initDumpStep = 500' and 'dataSetCount = 1000', it will use the 500th dump step & the 1000 after that.
 *
 *        If there are not sufficient dump steps to achieve the requested max (say only 750 in above example) as many as possible are used.  Normalization accounts for this automatically.
 *
 *          Useful if some simulations were not given sufficient walltime and didn't reach "completion," just use 'dataSetCount' = max of that throughout the simulation set.
 *  */

//  2018.04.10 NB Notes:  Usage Guidelines
/*
 * The executable performs various quantity of interest (QoI) calculations based on input flags.  Two types of flags exist:
 *
 *      Computation Flag:   Denoting which QoI(s) to compute.
 *
 *              Current options: 'A' for all, 'G' for RDF, 'M' for MSD, & 'T' for select dendrimer time-series.
 *
 *                  Note:  MSD is relative to the initial dump step, 'initDumpStep', specified.
 *
 *      Type Flag:          Denoting for which type of particles to compute the specified QoI(s).
 *
 *              Current options: 'A' for all, '1' for type-1 only, '2' for type-2 only.
 *
 * */

//  2018.04.10 NB Notes:  Common Usage Pitfalls
/*
 *
 * Must specify correct import file name, typically 'dump.melt', or you generate a blank instantaneous files & C++ chokes trying to import them.
 *
 * */

//  2018.04.10 NB Notes:  Necessary Code Improvements / Bugs
/*
 *  ### Desired Features ###
 *
 * Interval specification would be useful.  At present, 'dataSetCount' determines how many dump steps to use directly after 'initDumpStep'.
 *
 * It would - in principle - be wiser to import and do computations directly from the larger dump file, rather than smaller instantaenous dump files.
 *
 *   This would require adjusting how positions are updated, currently handled by sorting by particle index as the instantaneous files are generated.
 *
 * */

#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<assert.h>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>
#include "particle.h"
#include "bincontainer.h"
#include "functions.h"
#include<sstream>
#include<string>
#include<stdlib.h>
#include <cstring>

using namespace std;
using namespace boost::program_options;

// Fundamental quantities, units:
double pi = 3.141592;//65358979323846264338327950288419716939937510582097494459230781640628620899862803482; // Pi.
double Na = 6.022141e+23;                                                                               // Avogadro's.
double unitlength = 56.0e-9;//0.34e-09;      // Unit length (radius of virus);

double CondensedCount, BridgedCount;
vector<PARTICLE> Particle2_List_Condensed, Particle2_List_Bridged;
vector<double> CondensedCount_VsTime_List, BridgedCount_VsTime_List; // List of the number of selected particles of type 2 vs. used dump step.

// Declaration of RDF/MSD functions:
void assess_Select_Particles(vector<PARTICLE> &, vector<PARTICLE> &, double, double, double);

void
compute_gr_11(int, vector<BINCONTAINER> &, unsigned int, vector<PARTICLE> &, vector<PARTICLE> &, double, double, double,
              double, double, double);

void
compute_gr_12(int, vector<BINCONTAINER> &, unsigned int, vector<PARTICLE> &, vector<PARTICLE> &, double, double, double,
              double, double, double);

void
compute_gr_select_22(int, vector<BINCONTAINER> &, unsigned int, vector<PARTICLE> &, vector<PARTICLE> &, double, double,
                     double, double, double, double, int);

void compute_msd(vector<double> &, vector<PARTICLE> &, double, double, double);

int main(int argc, const char *argv[]) {

    /*Replacable variables*/
    string ligandNumText = "USERINPUT_LIGAND_NUMBER";
    string npChargeText = "USERINPUT_NP_CHARGE";
    string ligandChargeText = "DEVINPUT_LIGAND_CHARGE";
    string saltText = "USERINPUTSALT";
    string mpiText = "NODESIZE";

    // big for NP, small for ligand
    int n, Q;       //user input
    double c;       //user input
    double D, d;    //dev input, but could be made user input
    int q;          //dev input

    //  Declaring Type 1 particle {count, density, & diameter (forwarded but not used)}:
    long double Particle1_RealDensity, Particle2_RealDensity;
    unsigned int Particle1_Count, Particle2_Count;
    //  Input the molar concentrations (densities):
    Particle1_Count = 108;
    Particle1_RealDensity = .000000037 * 20;
    long double Particle1_Density = Na * 1000 * Particle1_RealDensity * pow(unitlength, 3);
    double Particle1_Diameter = 1;
//    double Particle2_Diameter = (6.7 / 56);
    //  Compute and report the resulting reduced box length, assign it to each component box length (for cubic):
    long double box_length = pow(Particle1_Count / (1000 * Na * Particle1_RealDensity), 1.0 / 3.0) /
                             unitlength; // 1000 because Molar is per dm^3.
    cout << "Computed box length (unitless): " << box_length << endl;
    cout << "\tComputed number of Type 1 particles in the box: " << Particle1_Count << endl;
    double bx = box_length, by = box_length, bz = box_length;
    //  Input the bin width intended:
    double bin_width = 0.005;//0.00005;

    //  Flag denoting computation type ('A' = all, '1' = 1-1 only, '2' = 1-2 & select 2-2 only).
    char computationFlag, typeFlag;
    //  Details on the datasets to be used (first dumpstep, number of directly subsequent steps):
    int initDumpStep, dataSetCount;
    //  The fold higher count of Type 2 particles:
    int stoichiometry;

    // Specify variables via command line (-X x):
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("Qnp,Q", value<int>(&Q)->default_value(-1500), "Q in e")
            ("NLigand,n", value<int>(&n)->default_value(25))
            ("Salt,c", value<double>(&c)->default_value(0.150), "c in Molars")
            ("qnp,q", value<int>(&q)->default_value(35), "q in e")
            ("NPDiameter,D", value<double>(&D)->default_value(56), "D in nm")
            ("LDiameter,d", value<double>(&d)->default_value(6.7), "d in nm")
            ("stoichiometry,x", boost::program_options::value<int>(&stoichiometry)->default_value(25),
             "Define which computation should be performed, 'A', 'G', 'M', 'T' only.")
            ("whichQuantity,q", boost::program_options::value<char>(&computationFlag)->default_value('G'),
             "Specify which computation should be performed, 'A', 'G' (RDF), 'M' (MSD), 'T' (Time-series) only.")
            ("whichTypes,t", boost::program_options::value<char>(&typeFlag)->default_value('1'),
             "Specify which type-type computations should be performed, 'A', '1', or '2' only.")
            ("initDumpStep,i", boost::program_options::value<int>(&initDumpStep)->default_value(0),
             "Specify the initial dump step to be used (dump step, not timestep).")
            ("dataSetCount,N", boost::program_options::value<int>(&dataSetCount)->default_value(100),
             "Specify the number of subsequent datasets to use after the initial dump step.");


    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);


    /*************************Post processing*************************/

    //  Defining the type 2 count and concentrations relative to the stoichiometry (w.r.t. type 1):
    stoichiometry = 0;
    Particle2_Count = stoichiometry * Particle1_Count;
    Particle2_RealDensity = stoichiometry * Particle1_RealDensity;
    cout << "\tComputed number of Type 2 particles in the box: " << Particle2_Count << endl;
    long double Particle2_Density = Na * 1000 * Particle2_RealDensity * pow(unitlength, 3);

    cout << "\nPreliminary quantities provided & computed.  Program begins.\n\n";

    //  Flags denoting if necessary to generate instantaneous dump files from larger movie file (and to delete afterwards):
    char createInstFilesFlag = 'y';
    char deleteInstFilesFlag = 'n';

    // Remove the previous and create the new instantaneous dump files (from the master movie 'dump.melt'):
    if (createInstFilesFlag == 'y') {
        if (boost::filesystem::remove_all("dumpfiles") == 0)
            cout << "\n\tError deleting instantaneous dump files directory ";
        else cout << "Pre-existing instantaneous dump files folder deleted successfully." << endl;
        //  Create the directory that will store instantaneous dump files & populate it:
        boost::filesystem::create_directory("dumpfiles");
        generate_inst_dump_files(Particle1_Count, Particle2_Count, initDumpStep, dataSetCount);
    }

    //  Declare and initialize the g(r) histogram list.
    vector<PARTICLE> dummy_particle_list, Particle1_List, Particle2_List;
    vector<BINCONTAINER> gr11, gr12, gr22;
    vector<double> MSD1, MSD2, MSD2_Select;
    if (computationFlag == 'A' || computationFlag == 'G') // [RDF-only]
    {
        if (typeFlag == 'A') {
            compute_gr_11(0, gr11, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, bin_width,
                          Particle1_RealDensity, Particle2_RealDensity);
            compute_gr_12(0, gr12, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, 0.1 * bin_width,
                          Particle1_RealDensity, Particle2_RealDensity);
            compute_gr_select_22(0, gr22, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, 0.1 * bin_width,
                                 Particle1_RealDensity, Particle2_RealDensity, initDumpStep);
        } else if (typeFlag == '1') {
            compute_gr_11(0, gr11, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, bin_width,
                          Particle1_RealDensity, Particle2_RealDensity);
        } else if (typeFlag == '2') {
            compute_gr_select_22(0, gr22, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, 0.1 * bin_width,
                                 Particle1_RealDensity, Particle2_RealDensity, initDumpStep);
        }
        cout << "\n RDF g(r) computation initialized, bins created.\n" << endl;
    }

    int actualDataSetCount = 0;  // NB added for normalization, in case a sim doesn't produce the requested max # of datasets (i.e. if it didn't finish).

    //  Import the data from each dataset and use it:
    for (int i = 0;
         i < dataSetCount; i++) { //  Read in the coordinates from instantaneous dump files, skipping header lines.
        //  Note, this 'instStream' is input, not the output "instStream" in 'generate_inst_dump_files(~)'.
        vector<VECTOR3D> newParticle1_Positions, newParticle2_Positions; // New particle positions, temporary for updating.
        int col1, col2;
        double col3, col4, col5;
        // Read data from a single specified file (within a for-loop iterating over all files).

        std::string fileName = "dumpfiles/";
        int fileNumber = initDumpStep + i;
        fileName = fileName + std::to_string(fileNumber) + ".melt";
        // If using 'split' bash (leading zeroes in timestep), verify %04 matches digit count.
        ifstream instStream(fileName, ios::in);
        if (!instStream)    // Verify the file could be opened.
        {
            cout << "dumpfiles" << fileNumber << ".melt could not be opened." << endl;
            continue;
        } else                // Warn the file could not be opened.
        {
            cout << "Opened file dumpfiles" << fileNumber << ".melt successfully." << endl;
            actualDataSetCount++;
        }

        // Skipping the first 9 lines in a standard LAMMPS output describing aspects of the simulation and step.
        for (int j = 1; j <= 9; j++) {
            string dummyline;
            getline(instStream, dummyline);
        }

        // Import the first (for constructing particle lists(s)) or latter (for updating particle properties) data files:
        while (instStream >> col1 >> col2 >> col3 >> col4 >> col5) {
            // Particle masses, charge (and diameters) do not matter, former set to zero.
            PARTICLE myparticle = PARTICLE(col1, Particle1_Diameter, 0, 0, VECTOR3D(col3, col4, col5),
                                           VECTOR3D(col3, col4, col5), bx, by, bz);
            //  If it's the initial data file, construct all particles and ascribe both position & initial position (same).
            if (i == 0) {
                if (col2 == 1) Particle1_List.push_back(myparticle);
                if (col2 == 2) Particle2_List.push_back(myparticle);
            } else { //  Else, construct a temporary list of the position vectors to simply update existing particles' positions:
                if (col2 == 1) newParticle1_Positions.push_back(VECTOR3D(col3, col4, col5));
                if (col2 == 2) newParticle2_Positions.push_back(VECTOR3D(col3, col4, col5));
            }
        }

        // Verify the imported information contained an initial or updated position for all particles as expected:
        if (i == 0) {
            //  Report and ensure the correct number of each have been constructed:
            cout << "Constructed number of particles of Type 1: " << Particle1_List.size() << endl;
            cout << "Constructed number of particles of Type 2: " << Particle2_List.size() << endl;
            assert(Particle1_Count ==
                   Particle1_List.size());   // Verify # of particles of each type imported are as expected.
            assert(Particle2_Count == Particle2_List.size());
        } else {
            //  Report and ensure the correct number of each have been constructed:
            cout << "Imported updates for particles of Type 1: " << newParticle1_Positions.size() << endl;
            cout << "Imported updates for particles of Type 2: " << newParticle2_Positions.size() << endl;
            assert(Particle1_Count ==
                   newParticle1_Positions.size());   // Verify all new positions were imported as expected.
            assert(Particle2_Count == newParticle2_Positions.size());
        }

        // If this is not the first data file, update particles' positions (using temp variable 'newParticleI_Positions'):
        if (i >= 1) { // This works because all instantaneous dump files have been sorted by particle index.
            for (unsigned int k = 0; k < Particle1_List.size(); k++) {
                Particle1_List[k].posvec = newParticle1_Positions[k];
            }
            for (unsigned int k = 0; k < Particle2_List.size(); k++) {
                Particle2_List[k].posvec = newParticle2_Positions[k];
            }
        }

        // Use the data (populate bins for RDF, compute MSD) for this dump file (pending requested QoI, types):
        if (computationFlag == 'A') {
            if (typeFlag == 'A') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                compute_gr_11(1, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                              Particle1_Density, 0 * Particle2_Density);
                compute_gr_12(1, gr12, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, 0.1 * bin_width,
                              Particle1_Density, Particle2_Density);
                compute_gr_select_22(1, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                     0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
                //compute_msd(MSD1, Particle1_List, bx, by, bz);
                //compute_msd(MSD2, Particle2_List, bx, by, bz);
                //compute_msd(MSD2_Select, Particle2_List_Condensed, bx, by, bz);
            } else if (typeFlag == '1') {
                compute_gr_11(1, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                              Particle1_Density, 0 * Particle2_Density);
                //compute_msd(MSD1, Particle1_List, bx, by, bz);
            } else if (typeFlag == '2') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                compute_gr_select_22(1, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                     0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
                //compute_msd(MSD2, Particle2_List, bx, by, bz);
                //compute_msd(MSD2_Select, Particle2_List_Condensed, bx, by, bz);
            }
        } else if (computationFlag == 'G') {
            if (typeFlag == 'A') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                compute_gr_11(1, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                              Particle1_Density, 0 * Particle2_Density);
                compute_gr_12(1, gr12, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, 0.1 * bin_width,
                              Particle1_Density, Particle2_Density);
                compute_gr_select_22(1, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                     0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
            } else if (typeFlag == '1') {
                compute_gr_11(1, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                              Particle1_Density, 0 * Particle2_Density);
            } else if (typeFlag == '2') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                compute_gr_select_22(1, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                     0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
            }
        } else if (computationFlag == 'M') {
            if (typeFlag == 'A') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                //compute_msd(MSD1, Particle1_List, bx, by, bz);
                //compute_msd(MSD2, Particle2_List, bx, by, bz);
                //compute_msd(MSD2_Select, Particle2_List_Condensed, bx, by, bz);
            } else if (typeFlag == '1') {
                //compute_msd(MSD1, Particle1_List, bx, by, bz);
            } else if (typeFlag == '2') {
                assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
                //compute_msd(MSD2, Particle2_List, bx, by, bz);
                //compute_msd(MSD2_Select, Particle2_List_Condensed, bx, by, bz);
            }
        } else if (computationFlag == 'T') assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);
    }

    // Report the number of data sets actually used vs. the number requested:
    if (actualDataSetCount == dataSetCount)
        cout << "All (" << dataSetCount << ") requested datasets imported successfully." << endl;
    else
        cout << "Warning:  only " << actualDataSetCount << " datasets were imported out of the total " << dataSetCount
             << " requested." << endl;

    // Delete the intermediate, instantaneous dump files (if requested):
    if (deleteInstFilesFlag == 'y') {
        boost::filesystem::remove_all("dumpfiles");
        cout << "\tDeleting instantaneous dump steps directory now that computation is complete." << endl;
    }

    //Particle1_List.resize(Particle1_Count);
    //Particle2_List.resize(Particle2_Count);
    // Normalize (and internally output) the RDF computation results:
    if (computationFlag == 'A' || computationFlag == 'G') {
        if (typeFlag == 'A') {
            compute_gr_11(2, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                          Particle1_Density, 0 * Particle2_Density);
            compute_gr_12(2, gr12, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, 0.1 * bin_width,
                          Particle1_Density, Particle2_Density);
            compute_gr_select_22(2, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                 0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
        } else if (typeFlag == '1') {
            compute_gr_11(2, gr11, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width,
                          Particle1_Density, 0 * Particle2_Density);
        } else if (typeFlag == '2') {
            compute_gr_select_22(2, gr22, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz,
                                 0.1 * bin_width, Particle2_Density, Particle2_Density, initDumpStep);
        }
    }

/*  // Output the MSD(t) data files (handled internally already for RDF):
  if (computationFlag == 'A' || computationFlag == 'M')
  {
    if (typeFlag == 'A')
    {
      ofstream msdStream1("MSD_V.out", ios::out);
      for (unsigned int i = 0; i < MSD1.size(); i++)  msdStream1 << (initDumpStep + i) << "\t" << MSD1[i] << endl;
      msdStream1.close();
      ofstream msdStream2("MSD_D.out", ios::out);
      for (unsigned int i = 0; i < MSD2.size(); i++)  msdStream2 << (initDumpStep + i) << "\t" << MSD2[i] << endl;
      msdStream2.close();
      ofstream msdStream2_Condensed("MSD_D_Condensed.out", ios::out);
      for (unsigned int i = 0; i < MSD2_Select.size(); i++)  msdStream2_Condensed << (initDumpStep + i) << "\t" << MSD2_Select[i] << endl;
      msdStream2_Condensed.close();
    }
    if (typeFlag == '1')
    {
      ofstream msdStream1("MSD_V.out", ios::out);
      for (unsigned int i = 0; i < MSD1.size(); i++)  msdStream1 << (initDumpStep + i) << "\t" << MSD1[i] << endl;
      msdStream1.close();
    }
    if (typeFlag == '2')
    {
      ofstream msdStream2("MSD_D.out", ios::out);
      for (unsigned int i = 0; i < MSD2.size(); i++)  msdStream2 << (initDumpStep + i) << "\t" << MSD2[i] << endl;
      msdStream2.close();
      ofstream msdStream2_Condensed("MSD_D_Condensed.out", ios::out);
      for (unsigned int i = 0; i < MSD2_Select.size(); i++)  msdStream2_Condensed << (initDumpStep + i) << "\t" << MSD2_Select[i] << endl;
      msdStream2_Condensed.close();
    }
  }*/

    // Output the select dendrimer time-series data files:
    if ((computationFlag == 'A' || computationFlag == 'T') && (typeFlag == 'A' || typeFlag == '2')) {
        ofstream SelectCountStream("SelectDendrimer_Count_VsTime_dr=1.05.out", ios::out);
        for (unsigned int i = 0; i < CondensedCount_VsTime_List.size(); i++)
            SelectCountStream << initDumpStep + i << "\t" << CondensedCount_VsTime_List[i] << "\t"
                              << BridgedCount_VsTime_List[i] << endl;
        // Abscissae are by the used dump steps, i.e. if you start at step 200, "1" is 200.
        SelectCountStream.close();
        cout << "Select dendrimer time-series file output complete." << endl;
    }

    cout << "Program complete." << endl;
    return 0;
}
// End of main
