
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options.hpp>


using namespace std;
using namespace boost::program_options;



int main(int argc, const char *argv[]) {

    /*User Inputs*/
    //string ligandNumText = "USERINPUT_LIGAND_NUMBER";
    //string npChargeText = "USERINPUT_NP_CHARGE";
    //string mpiText = "NODESIZE";

    string ligandNumFileNameText = "USERINPUT_LIGAND_FILENAME";
    string saltText = "USERINPUTSALT";
    string dataSetCountText = "USERINPUT_DATASETCOUNT";

    // Dev inputs
    string vlp1ChargeText = "DEVINPUT_VLP1_CHARGE";
    string vlp2ChargeText = "DEVINPUT_VLP2_CHARGE";
    string ligandChargeText = "DEVINPUT_LIGAND_CHARGE";

    const double pi = 3.141593;
    // big for NP, small for ligand
    int n, Q_vlp1, Q_vlp2;       //dev input
    int e2k2ratio; //user input
    double c;       //user input
    double D, d;    //dev input, but could be made user input
    int q;          //dev input

    int mpiProcs;

    //  Details on the datasets to be used (first dumpstep, number of directly subsequent steps):
    int initDumpStep, dataSetCount;


    // Specify variables via command line (-X x):
    options_description desc("Usage:\nrandom_mesh <options>");
    desc.add_options()
            ("help,h", "print usage message")
            ("Qvlp1,E", value<int>(&Q_vlp1)->default_value(-1500), "Q_E2 in e")
            ("Qvlp2,K", value<int>(&Q_vlp2)->default_value(-600), "Q_K2 in e")
            ("e2k2ratio,r", value<int>(&e2k2ratio)->default_value(1), "e2:k2 in multiples")
            ("NLigand,n", value<int>(&n)->default_value(100))
            ("Salt,c", value<double>(&c)->default_value(0.150), "c in Molars")
            ("qnp,q", value<int>(&q)->default_value(45), "q in e")
            ("NPDiameter,D", value<double>(&D)->default_value(56), "D in nm")
            ("LDiameter,d", value<double>(&d)->default_value(6.7), "d in nm")
            ("Mpi_procs,m", value<int>(&mpiProcs)->default_value(1), "Number of MPI procs for Lammps")
            ("initDumpStep,i", boost::program_options::value<int>(&initDumpStep)->default_value(0),
             "Specify the initial dump step to be used (dump step, not timestep).")
            ("dataSetCount,N", boost::program_options::value<int>(&dataSetCount)->default_value(150),
             "Specify the number of subsequent datasets to use after the initial dump step.");
    //hard code the ratio

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);


    /*************************Preprocessing*************************/
    /*Open a new the template file*/


    /******************computations of dev variables****************/

    double devInputQvlp1Charge = Q_vlp1 * (exp((1.64399 * 56) / sqrt(1 / c)) / (1 + ((1.64399 * 56) / sqrt(1 / c))));
    devInputQvlp1Charge = devInputQvlp1Charge * (1.6018 * 1e-19) / (1.60074 * 1e-19);
    
    double devInputQvlp2Charge = Q_vlp2 * (exp((1.64399 * 56) / sqrt(1 / c)) / (1 + ((1.64399 * 56) / sqrt(1 / c))));
    devInputQvlp2Charge = devInputQvlp2Charge * (1.6018 * 1e-19) / (1.60074 * 1e-19);
    
    double devInputLigandCharge = q * (exp((1.64399 * 6.7) / sqrt(1 / c)) / (1 + ((1.64399 * 6.7) / sqrt(1 / c))));
    devInputLigandCharge = devInputLigandCharge * (1.6018 * 1e-19) / (1.60074 * 1e-19);

    double userInputSalt = 1 / (1e3 * sqrt(1 / (2 * 4 * pi * 0.714295 * (6.022 * 1e23) * c)));
    //double userInputSalt = 3.2879795708993616*1e9/sqrt(1/c);
    userInputSalt = userInputSalt * 56 * 1e-9;
    
    
    std::string initcoords = "initCoords_" + std::to_string(n);
    if(e2k2ratio==1)
        initcoords = initcoords+"x_1-1";
    else if(e2k2ratio == 4)
        initcoords = initcoords+"x_1-4";
    initcoords = initcoords + ".assembly";
    //cout<<"\n ===>"<<initcoords;
    ofstream inputScript("in.lammps", ios::trunc);
    if (inputScript.is_open()) {

        /*Open the template file*/
        string line;
        ifstream inputTemplate("infiles/in.lammps.template", ios::in);
        if (inputTemplate.is_open()) {
            while (getline(inputTemplate, line)) {
//                std::size_t found = line.find(ligandNumText);
//                if (found != std::string::npos)
//                    line.replace(found, ligandNumText.length(),  std::to_string(n));

//                found = line.find(npChargeText);
//                if (found != std::string::npos)
//                    line.replace(found, npChargeText.length(), std::to_string(userInputNpCharge));

                // Added for bnsl
                std::size_t found = line.find(vlp1ChargeText);
                if (found != std::string::npos)
                    line.replace(found, vlp1ChargeText.length(), std::to_string(devInputQvlp1Charge));
                found = line.find(vlp2ChargeText);
                if (found != std::string::npos)
                    line.replace(found, vlp2ChargeText.length(), std::to_string(devInputQvlp2Charge));

                found = line.find(ligandChargeText);
                if (found != std::string::npos)
                    line.replace(found, ligandChargeText.length(), std::to_string(devInputLigandCharge));

                found = line.find(saltText);
                if (found != std::string::npos)
                    line.replace(found, saltText.length(), std::to_string(userInputSalt));

                //unsure
                found = line.find(dataSetCountText);
                if (found != std::string::npos)
                    line.replace(found, dataSetCountText.length(), std::to_string(dataSetCount));

               found = line.find(ligandNumFileNameText);
               if (found != std::string::npos)
                   line.replace(found, ligandNumFileNameText.length(), initcoords);

                inputScript << line << endl;
            }
            inputTemplate.close();
        } else cout << "Unable to open the template input script" << endl;
        inputScript.close();
    } else cout << "Unable create a input Script" << endl;


    /*************************Lammps Call*************************/
/*
    //string lammpsExeCMD = "aprun -n NODESIZE lmp_mpi < in.lammps"; //This is Bigred and New Lammps
    //string lammpsExeCMD = "mpirun -n NODESIZE lmp_g++ < in.lammps"; //This is RedHat
    string lammpsExeCMD = "mpiexec -n NODESIZE lmp_mpi < in.lammps"; //This is windows and newest Lammps
    //string lammpsExeCMD = "aprun -n NODESIZE lmp_xe6  < in.lammps"; //This is Bigred and OLD lammps
    std::size_t found = lammpsExeCMD.find(mpiText);
    if (found != std::string::npos)
        lammpsExeCMD.replace(found, mpiText.length(), std::to_string(mpiProcs));

    //Creating the char array
    char lammpsExeCMD_array[lammpsExeCMD.length() + 1];
    strcpy(lammpsExeCMD_array, lammpsExeCMD.c_str());

    //System call
    //system("ls");
    int returnedCode;
    //Checking if processor is available
    if (!system(NULL))
        exit(EXIT_FAILURE);

    returnedCode = system(lammpsExeCMD_array);
    cout << "The value returned was: " << returnedCode << endl;

    if(returnedCode!=0)
        exit(EXIT_FAILURE);
*/
    cout << "Preprocessing completed." << endl;
    return 0;
}
// End of main
