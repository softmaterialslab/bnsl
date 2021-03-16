/* 
 * Feb, 2021
 * get PCF for type 1 particles
 * rewrite, streamline the code and organize content for future updates
 *
 * need to change double to long double: unitlength, bx by changed type
*/

#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<assert.h>
#include<boost/filesystem/operations.hpp>
#include<boost/program_options.hpp>
#include"particle.h"
#include"bincontainer.h"
#include"functions.h"
#include<sstream>
#include<string>
#include<stdlib.h>
#include<cstring>

using namespace std;
using namespace boost::program_options;

// Fundamental quantities, units:
double pi = 3.141592;//65358979323846264338327950288419716939937510582097494459230781640628620899862803482; // Pi.
double Na = 6.022141e+23;	// Avogadro's.
long double unitlength; 	// Unit length (diameter of virus in meters);

// Declaration of RDF functions:
void compute_self_gr(int, vector<BINCONTAINER>&, unsigned int, vector<PARTICLE>&, long double, long double, long double, double, long double, ofstream&);
void compute_gr_12(int, vector<BINCONTAINER> &, unsigned int, vector<PARTICLE> &, vector<PARTICLE> &, long double, long double, long double, double, long double, long double);

/* Later dev:
void assess_Select_Particles(vector<PARTICLE>&, vector<PARTICLE>&, double, double, double);
double CondensedCount, BridgedCount;
vector<PARTICLE> Particle2_List_Condensed, Particle2_List_Bridged;
vector<double> CondensedCount_VsTime_List, BridgedCount_VsTime_List; // List of the number of selected particles of type 2 vs. used dump step.
*/

int main(int argc, const char *argv[]) {
    // User defined Inputs
    int stoichiometry;      // The fold higher count of Type 2 particles

    unsigned int Particle1_Count, Particle2_Count;
    long double Particle1_RealDensity, Particle2_RealDensity;
    long double D, d;

    //  Flag denoting computation type ('A' = all, '1' = 1-1 only, '2' = 1-2 & select 2-2 only).
    char computationFlag, typeFlag;
    //  Details on the datasets to be used (first dumpstep, number of directly subsequent steps)
    int initDumpStep, dataSetCount;
    
    // Specify variables via command line (-X -x)
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("VLP1count,e", value<unsigned int>(&Particle1_Count)->default_value(50), "Number of particles of type 1")
            ("VLP1diameter,D", value<long double>(&D)->default_value(56), "Diameter of type 1 particles (in nm)")
            ("VLP2diameter,d", value<long double>(&d)->default_value(56), "d in nm")
            ("stoichiometry,x", boost::program_options::value<int>(&stoichiometry)->default_value(1), "Relative ratio of two types (check)")
            ("whichComp,q", boost::program_options::value<char>(&computationFlag)->default_value('G'), "Specify computation: 'G' (RDF), 'T' (Time-series)")
            ("whichTypes,t", boost::program_options::value<char>(&typeFlag)->default_value('1'), "Specify type-type: 'A', '1', or '2'")
            ("initDumpStep,i", boost::program_options::value<int>(&initDumpStep)->default_value(0), "Specify the initial dump step (not timestep)")
            ("dataSetCount,N", boost::program_options::value<int>(&dataSetCount)->default_value(2001), "Specify the number of subsequent time samples");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    unitlength = D/1000000000;
    double Particle1_Diameter = 1;			// in reduced units of particle 1 diameter (not used)

    //  Input the molar concentrations (densities)
    Particle1_RealDensity = (370e-9) / (1+stoichiometry);		// 370 nM is total density; half of this is particle 1 density

    //  Compute density in reduced units (1L = 0.001 m^3)
    long double Particle1_Density = Particle1_RealDensity * Na * 1000 * pow(unitlength, 3);	

    //  Compute the resulting box length, assign it to cubic box dimensions:
    long double box_length = pow(Particle1_Count / Particle1_Density, 1.0 / 3.0); 
    long double bx = box_length, by = box_length, bz = box_length;

    cout << "Number of Type 1 particles in the box: " << Particle1_Count << endl;
    cout << "Computed box length (reduced units): " << box_length << endl;

    //  Input the bin width intended:
    double bin_width = 0.005;	// change to smaller if needed

    //  Defining the type 2 count and concentrations relative to the stoichiometry (w.r.t. type 1):
    Particle2_Count = stoichiometry * Particle1_Count;
    Particle2_RealDensity = stoichiometry * Particle1_RealDensity;
    cout << "Computed number of Type 2 particles in the box: " << Particle2_Count << endl;
    long double Particle2_Density = Particle2_RealDensity * Na * 1000 * pow(unitlength, 3);

    cout << "\nPreliminary quantities provided & computed. Postprocessing begins.\n";

    /*************************Post processing actions*************************/

    //  Flags denoting if necessary to generate instantaneous dump files from larger movie file (and to delete afterwards):
    char createInstFilesFlag = 'y';
    char deleteInstFilesFlag = 'n';

    //  Remove the previous and create the new instantaneous dump files (from the master movie 'dump.melt'):
    if (createInstFilesFlag == 'y') {
        if (boost::filesystem::remove_all("dumpfiles") == 0)
            cout << "\nNo files in dumpfiles directory\n";
        else cout << "Pre-existing instantaneous dumpfiles directory deleted\n" << endl;
        //  Create the directory that will store instantaneous dump files & populate it:
        boost::filesystem::create_directory("dumpfiles");
	    cout << "Generating samples..." << endl;
        generate_inst_dump_files(Particle1_Count, Particle2_Count, initDumpStep, dataSetCount);
    }

    //  Declare and initialize the g(r) histogram list.
    vector<PARTICLE> dummy_particle_list, Particle1_List, Particle2_List;
    vector<BINCONTAINER> gr11, gr12, gr22;
	 
    ofstream grStream11("outfiles/gr_EE_dr=0.005.out", ios::out);
    ofstream grStream22("outfiles/gr_KK_dr=0.005.out", ios::out);

    cout << "RDF g(r) computation initialized" << endl;

    if (computationFlag == 'G') // [RDF-only]
    {
        if (typeFlag == 'A') {
            compute_self_gr(0, gr11, 0, dummy_particle_list, bx, by, bz, bin_width, Particle1_RealDensity, grStream11);
            compute_self_gr(0, gr22, 0, dummy_particle_list, bx, by, bz, bin_width, Particle2_RealDensity, grStream22);
            //compute_gr_12(0, gr12, 0, dummy_particle_list, dummy_particle_list, bx, by, bz, bin_width, Particle1_Density, Particle2_Density);
        }
	else if (typeFlag == '1') {
            compute_self_gr(0, gr11, 0, dummy_particle_list, bx, by, bz, bin_width, Particle1_RealDensity, grStream11);
        }
	else if (typeFlag == '2') {
            compute_self_gr(0, gr22, 0, dummy_particle_list, bx, by, bz, bin_width, Particle2_RealDensity, grStream22);
        }
    }

    int actualDataSetCount = 0;  // in case a sim doesn't produce the requested max # of datasets (i.e. if it didn't finish).

    //  Import the data from each sample and use it
    //  Read data from a single specified file (within a for-loop iterating over all files).
    
    for (int i = 0;  i < dataSetCount; i++) { // Read in the coordinates from instantaneous dump files, skipping header lines.
        vector<VECTOR3D> newParticle1_Positions, newParticle2_Positions; // New particle positions, temporary for updating.
        int col1, col2;
        double col3, col4, col5;

        std::string fileName = "dumpfiles/";
        int fileNumber = initDumpStep + i;
        fileName = fileName + std::to_string(fileNumber) + ".melt"; // If using 'split' bash (leading zeroes in timestep), verify %04 matches digit count
        
        //  Note, this 'instStream' is input, not the output "instStream" in 'generate_inst_dump_files(~)'.
        ifstream instStream(fileName, ios::in);
        if (!instStream) {
            cout << "dumpfiles" << fileNumber << ".melt could not be opened." << endl;
            continue;
        } 
        
        else {
            if (fileNumber%100==0) 
					cout << "Opened sample " << fileNumber << ".melt successfully." << endl;
            actualDataSetCount++;
        }

        //  Skipping the first 9 lines in a standard LAMMPS output describing aspects of the simulation and step.
        for (int j = 1; j <= 9; j++) {
            string dummyline;
            getline(instStream, dummyline);
        }

        //  Import the first sample for constructing particle list; latter samples update particle properties
        while (instStream >> col1 >> col2 >> col3 >> col4 >> col5) {
            //  Particle mass, charge set to zero (does not matter for PCF)
            PARTICLE myparticle = PARTICLE(col1, Particle1_Diameter, 0, 0, VECTOR3D(col3, col4, col5), VECTOR3D(col3, col4, col5), bx, by, bz);
            //  If it's the initial data file, construct all particles and ascribe both position & initial position (same).
            if (i == 0) {
                if (col2 == 1) Particle1_List.push_back(myparticle);
                if (col2 == 2) Particle2_List.push_back(myparticle);
            } else { //  Else, construct a temporary list of the position vectors to simply update existing particles' positions:
                if (col2 == 1) newParticle1_Positions.push_back(VECTOR3D(col3, col4, col5));
                if (col2 == 2) newParticle2_Positions.push_back(VECTOR3D(col3, col4, col5));
            }
        }

        //  Verify the imported information 
        if (i == 0) {
            //  Report and ensure the correct number of each have been constructed:
            cout << "\tConstructed number of particles of Type 1: " << Particle1_List.size() << " Type 2: " << Particle2_List.size() << endl;
            assert(Particle1_Count == Particle1_List.size());   	// Verify # of particles of each type imported are as expected.
            assert(Particle2_Count == Particle2_List.size());
        } else if (fileNumber%100==0){
            //  Report and ensure the correct number of each have been constructed:
            cout << "\tImported updates for particles of Type 1: " << newParticle1_Positions.size() << " Type 2: " << newParticle2_Positions.size() << endl;
            assert(Particle1_Count == newParticle1_Positions.size());   // Verify all new positions were imported as expected.
            assert(Particle2_Count == newParticle2_Positions.size());
        }

        //  If this is not the first data file, update particles' positions using temp variable 'newParticle1_Positions'
        if (i >= 1) { //  This works because all instantaneous dump files have been sorted by particle index.
            for (unsigned int k = 0; k < Particle1_List.size(); k++) {
                Particle1_List[k].posvec = newParticle1_Positions[k];
            }
            for (unsigned int k = 0; k < Particle2_List.size(); k++) {
                Particle2_List[k].posvec = newParticle2_Positions[k];
            }
        }

        //  Use the data (populate bins for RDF) for this dump file (pending requested QoI, types):
        if (computationFlag == 'G') {
            if (typeFlag == 'A') {
                compute_self_gr(1, gr11, actualDataSetCount, Particle1_List, bx, by, bz, bin_width, Particle1_Density, grStream11);
                compute_self_gr(1, gr22, actualDataSetCount, Particle2_List, bx, by, bz, bin_width, Particle2_Density, grStream22);
                //compute_gr_12(1, gr12, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width, Particle1_Density, Particle2_Density);
            } 
	    else if (typeFlag == '1') {
                compute_self_gr(1, gr11, actualDataSetCount, Particle1_List, bx, by, bz, bin_width, Particle1_Density, grStream11);
            } 
	    else if (typeFlag == '2') {
                compute_self_gr(1, gr22, actualDataSetCount, Particle2_List, bx, by, bz, bin_width, Particle2_Density, grStream22);
            }
        } 

    } // end of for loop

    // if (computationFlag == 'T') assess_Select_Particles(Particle1_List, Particle2_List, bx, by, bz);

    //  Report the number of data sets actually used vs the number requested
    if (actualDataSetCount == dataSetCount)
        cout << "All (" << dataSetCount << ") requested datasets imported successfully" << endl;
    else
        cout << "Warning:  only " << actualDataSetCount << " datasets were imported out of the total " << dataSetCount << " requested" << endl;

    // Delete the intermediate, instantaneous dump files (if requested)
    if (deleteInstFilesFlag == 'y') {
        cout << "\tDeleting instantaneous dump steps directory...";
        boost::filesystem::remove_all("dumpfiles");
        cout << "done" << endl;
    }

    //  Normalize output the RDF computation results:
    if (computationFlag == 'G') {
        if (typeFlag == 'A') {
            compute_self_gr(2, gr11, actualDataSetCount, Particle1_List, bx, by, bz, bin_width, Particle1_Density, grStream11);
            compute_self_gr(2, gr22, actualDataSetCount, Particle2_List, bx, by, bz, bin_width, Particle2_Density, grStream22);
            //compute_gr_12(2, gr12, actualDataSetCount, Particle1_List, Particle2_List, bx, by, bz, bin_width, Particle1_Density, Particle2_Density);
        } 
	else if (typeFlag == '1') {
            compute_self_gr(2, gr11, actualDataSetCount, Particle1_List, bx, by, bz, bin_width, Particle1_Density, grStream11);
        } 
	else if (typeFlag == '2') {
            compute_self_gr(2, gr22, actualDataSetCount, Particle2_List, bx, by, bz, bin_width, Particle2_Density, grStream22);
        }
    }

    return 0;
}
