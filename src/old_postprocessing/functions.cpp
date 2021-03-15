// This file contains the calculation steps expecting to be fed as arguments the binning apparatus (initially empty g(r)), the type-segregated coordinate sets, and the corresponding types' (Molar) densities.

#include "functions.h"
#include "bincontainer.h"
#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>      // For the 'sort()' function used in generating type-specific instantaneous dump files.

extern vector<PARTICLE> Particle2_List_Condensed, Particle2_List_Bridged;
extern double pi, Na, unitlength, CondensedCount, BridgedCount;

// Quantities used within functions, but that can't be re-declared or initialized each time:
extern vector<double> CondensedCount_VsTime_List, BridgedCount_VsTime_List; // List of the number of selected particles of type 2 vs. used dump step.

//  Overload the output operator for use with the VECTOR3D class:
ostream& operator<<(ostream& os, VECTOR3D vec)
{
  os << vec.x << setw(15) << vec.y << setw(15) << vec.z;
  return os;
}

// A comparator function to sort particles by index (thus also sorting by type, if correlated):
bool is_not_digit(char c)
{
    return !std::isdigit(c);
}
bool numeric_string_compare(const std::string& s1, const std::string& s2)
{
    std::string::const_iterator it1 = s1.begin(), it2 = s2.begin();

    if (std::isdigit(s1[0]) && std::isdigit(s2[0])) {
        int n1, n2;
        std::stringstream ss(s1);
        ss >> n1;
        ss.clear();
        ss.str(s2);
        ss >> n2;

        if (n1 != n2) return n1 < n2;

        it1 = std::find_if(s1.begin(), s1.end(), is_not_digit);
        it2 = std::find_if(s2.begin(), s2.end(), is_not_digit);
    }

    return std::lexicographical_compare(it1, s1.end(), it2, s2.end());
}

//  Generate instantaneous (single timestep) files from the larger movie (dump) file:
void generate_inst_dump_files(int Particle1_Count, int Particle2_Count, int initDumpStep, int dataSetCount)
{
  int particleCount = Particle1_Count + Particle2_Count;
  vector<string> instDumpLines;             //  List of all lines to be included in the instantaneous dump file.
  ifstream movieStream;                     //  Master dumpfile from which to extract instantaneous dump files.
  movieStream.open("outfiles/dump.melt");            //  Opens the master dumpfile.
  //char instFileName[100];                   //  Buffers to be filled with the instantaneous file names.
    std::string instFileName = "dumpfiles/";
  ofstream instStream(instFileName, ios::out);  //  The stream that will be output to the all-types file.
  string tempString;                        //  A temporary string that will be updated as the file is read.
  int j = 0;                                //  Dummy (line/stream position) variable.

  double fileNumber;                        //  The file number (if rounded, dump step) being exported.
  //  While the movie stream hasn't reached its end, continue importing and exporting requested content:
  while (!movieStream.eof())
  {
    j++;                                                      //  Iterate to track (line) position of the stream.
    fileNumber = ((double)j/(9.0+particleCount) - 1.0);       //  The current dumpstep (net # lines / # per step).
      int fileNumberInt = (int) fileNumber;
    getline(movieStream, tempString);
    instDumpLines.push_back(tempString);
    if(((j%(9 + particleCount)) == 0) && (initDumpStep <= fileNumber) && (fileNumber < (initDumpStep + dataSetCount)))
    {


        instFileName = "dumpfiles/" + std::to_string(fileNumberInt) + ".melt";

        //sprintf(instFileName, "dumpfiles/%d.melt", int(fileNumber));  // Define file name.
      instStream.open(instFileName);                                                          // Open stream to file.
      sort(instDumpLines.end()-particleCount,instDumpLines.end(),numeric_string_compare);     // Sort by index.
      for(vector<string>::iterator it = instDumpLines.begin() ; it != instDumpLines.end(); ++it)
        instStream << *it << endl;
      instStream.close();
      cout << "\tGenerated index-sorted ensemble instantaneous dump file " << fileNumber << " from master movie." << endl;
    }
    if((j%(9 + particleCount)) == 0) instDumpLines.clear();
    if(fileNumber > (initDumpStep + dataSetCount)) break; // Breaks the loop if the fileNumber exceeds that requested.
  }
  movieStream.close();
}

// The following RDF functions differ minimally from one another.  E.g. the only change in "gr_22" is output file name.

// A function to assess the condensed and bridging type-2 particles (relative to type-1):
void assess_Select_Particles(vector<PARTICLE>& Particle1_List, vector<PARTICLE> &Particle2_List, double bx, double by, double bz)
{
    //  Iterate over all type 2 particles, selecting only those within a threshold distance from particle(s) of type 1:
    Particle2_List_Condensed.clear();   // Clear elements stored from the previous function call / dump step.
    Particle2_List_Bridged.clear();
    for (unsigned int i = 0; i < Particle2_List.size(); i++)        // Iterates over all type-2 (dendrimers).
    {
        int vCountNear = 0;
        for (unsigned int j = 0; j < Particle1_List.size(); j++)    //  Iterates over all type-1 (viruses) for each type-2 (dendrimer).
        {
            VECTOR3D r_vec = Particle2_List[i].posvec - Particle1_List[j].posvec;
            if (r_vec.x>bx/2) r_vec.x -= bx;
            if (r_vec.x<-bx/2) r_vec.x += bx;
            if (r_vec.y>by/2) r_vec.y -= by;
            if (r_vec.y<-by/2) r_vec.y += by;
            if (r_vec.z>bz/2) r_vec.z -= bz;
            if (r_vec.z<-bz/2) r_vec.z += bz;
            double r=r_vec.GetMagnitude();

            if(r < 1.05*(((56 + 6.7)/2)/56)) vCountNear++; // Record the number of type-1 particles a given type-2 particle is near.
        }
        if(vCountNear >= 1) Particle2_List_Condensed.push_back(Particle2_List[i]);
        if(vCountNear >= 2) Particle2_List_Bridged.push_back(Particle2_List[i]);
    }

    CondensedCount = Particle2_List_Condensed.size();
    BridgedCount = Particle2_List_Bridged.size();
    cout << "\t\tThe number of condensed dendrimer this timestep is: " << CondensedCount << endl;
    cout << "\t\tThe number of bridging dendrimer this timestep is: " << BridgedCount << endl;
    CondensedCount_VsTime_List.push_back(CondensedCount);
    BridgedCount_VsTime_List.push_back(BridgedCount);
}

//  Compute RDF g(r): pair correlation function for type 1-1 particles (say, virus-virus):
void compute_gr_11(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE>& Particle1_List, vector<PARTICLE>& Particle2_List, double bx, double by, double bz, double bin_width, double Particle1_Density, double Particle2_Density)
{
  if (stageFlag == 0)	// Initialize the ensemble g(r) bins.
  {
    ngr=0;		    // Number of datasets.
    int number_of_bins = int((bz/2.0)/bin_width);
    cout << "The number of bins is: " << number_of_bins << endl;
    gr.resize(number_of_bins);
    for (int i = 0; i < number_of_bins; i++)
    {
      gr[i].number = i+1;
      gr[i].width = bin_width;
      gr[i].position = (gr[i].width)*i;             // Assigns bin positions.  See note directly below as well.
      gr[i].population = 0.0;                       // NB moved moving median (center of bins) to post-normalization.
    }
  }
  if (stageFlag == 1)	// Bin all distances contained in the input file.
  { // To change between poly- and monodisperse, values denoted in right margin.
  for (unsigned int i = 0; i < Particle1_List.size(); i++)                      // For 1 type, change here (1).
    {
      for (unsigned int j = i+1; j < Particle1_List.size(); j++)                    // For 1 type, change here (2).
      {
        VECTOR3D r_vec = Particle1_List[i].posvec - Particle1_List[j].posvec;     // For 1 type, change here (3).
        if (r_vec.x>bx/2) r_vec.x -= bx;
        if (r_vec.x<-bx/2) r_vec.x += bx;
        if (r_vec.y>by/2) r_vec.y -= by;
        if (r_vec.y<-by/2) r_vec.y += by;
        if (r_vec.z>bz/2) r_vec.z -= bz;
        if (r_vec.z<-bz/2) r_vec.z += bz;
        double r=r_vec.GetMagnitude();
        if (r < bz/2.0 - 1)	// avoiding the g(r) calculation at the largest r
        {
          int bin_number = ceil((r/bin_width));
          gr[bin_number - 1].population = gr[bin_number - 1].population + 2;  // For 1 type, optionally change here (3.1).
        }
      }
    }
    cout << "\tDataset binning complete for RDF (type 1-1)." << endl;
  }
  if (stageFlag == 2)	// Normalize each bin count appropriately (see comments).
  {
    cout << endl << "RDF (type 1-1) calculation ends, beginning normalization & output." << endl;
    int number_of_bins = int((bz/2.0)/bin_width);
    ofstream grStream("outfiles/gr_VV_dr=0.005.out", ios::out);

    for (int b = 0; b < number_of_bins; b++)
    { //  Output moving median abscissae (rMedian) & normalize for g(r) using non-moving median shell volumes binned.
            double r = gr[b].position;
      double vol_bin = (4.0 / 3.0) * pi * ((pow(r + bin_width, 3) - pow(r, 3)));  // Concentric shells' interstitial volume.
      double nid = vol_bin * Particle1_Density;                       // For 1 type, change ParticleJ_Density (4).
      double rMedian = r + 0.5*bin_width;                             // Compute the moving median as final abscissae.
      gr[b].population = gr[b].population / Particle1_List.size();    // Normalize by the number of particles.
      gr[b].population = gr[b].population / ngr;                      // Normalize by the number of datasets.
      gr[b].population = gr[b].population / nid;                      // Normalize by the expected number in ideal gas.
      if (rMedian <= 4.3)                                             // Limiting to not show data at artificial cutoff
        grStream << rMedian << "\t" << gr[b].population << endl;
    }
    grStream.close();
    cout << "\tRDF (type 1-1) file output complete." << endl;
  }
}

//  Compute RDF g(r): pair correlation function for type 1-2 particles (say, virus-dendrimer):
void compute_gr_12(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE>& Particle1_List, vector<PARTICLE>& Particle2_List, double bx, double by, double bz, double bin_width, double Particle1_Density, double Particle2_Density)
{
  if (stageFlag == 0)	// Initialize the ensemble g(r) bins.
  {
    ngr=0;		    // Number of datasets.
    int number_of_bins = int((bz/2.0)/bin_width);
    cout << "The number of bins is: " << number_of_bins << endl;
    gr.resize(number_of_bins);
    for (int i = 0; i < number_of_bins; i++)
    {
      gr[i].number = i+1;
      gr[i].width = bin_width;
      gr[i].position = (gr[i].width)*i;             // Assigns bin positions.  See note directly below as well.
      gr[i].population = 0.0;                       // NB moved moving median (center of bins) to post-normalization.
    }
  }
  if (stageFlag == 1)	// Bin all distances contained in the input file.
  { // To change between poly- and monodisperse, values denoted in right margin.
    for (unsigned int i = 0; i < Particle1_List.size(); i++)                      // For 1 type, change here (1).
    {
      for (unsigned int j = 0; j < Particle2_List.size(); j++)                    // For 1 type, change here (2).
      {
        VECTOR3D r_vec = Particle1_List[i].posvec - Particle2_List[j].posvec;     // For 1 type, change here (3).
        if (r_vec.x>bx/2) r_vec.x -= bx;
        if (r_vec.x<-bx/2) r_vec.x += bx;
        if (r_vec.y>by/2) r_vec.y -= by;
        if (r_vec.y<-by/2) r_vec.y += by;
        if (r_vec.z>bz/2) r_vec.z -= bz;
        if (r_vec.z<-bz/2) r_vec.z += bz;
        double r=r_vec.GetMagnitude();
        if (r < bz/2.0)
        {
          int bin_number = ceil((r/bin_width));
          gr[bin_number - 1].population = gr[bin_number - 1].population + 1;  // For 1 type, optionally change here (3.1).
        }
      }
    }
    cout << "\tDataset binning complete for RDF (type 1-2)." << endl;
  }
  if (stageFlag == 2)	// Normalize each bin count appropriately (see comments).
  {
    cout << endl << "RDF (type 1-2) calculation ends, beginning normalization & output." << endl;
    int number_of_bins = int((bz/2.0)/bin_width);
    ofstream grStream("gr_VD_dr=0.0005.out", ios::out);

    for (int b = 0; b < number_of_bins; b++)
    { //  Output moving median abscissae (rMedian) & normalize for g(r) using non-moving median shell volumes binned.
      double r = gr[b].position;
      double vol_bin = (4.0 / 3.0) * pi * ((pow(r + bin_width, 3) - pow(r, 3)));  // Concentric shell interstitial volume.
      double nid = vol_bin * Particle2_Density;                       // For 1 type, change ParticleJ_Density (4).
      double rMedian = r + 0.5*bin_width;                             // Compute the moving median as final abscissae.
      gr[b].population = gr[b].population / Particle1_List.size();    // Normalize by the number of particles.
      gr[b].population = gr[b].population / ngr;                      // Normalize by the number of datasets.
      gr[b].population = gr[b].population / nid;                      // Normalize by the expected number in ideal gas.
      grStream << rMedian << "\t" << gr[b].population << endl;
      //grStream.close();
    }
    cout << "\tRDF (type 1-2) file output complete." << endl << endl;
  }
}

//  Compute RDF g(r): pair correlation function for select type 2-2 particles (say, condensed-only dendrimer-dendrimer):
void compute_gr_select_22(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE> &Particle1_List, vector<PARTICLE> &Particle2_List, double bx, double by, double bz, double bin_width, double Particle1_Density, double Particle2_Density, int initDumpStep)
{
  if (stageFlag == 0)	// Initialize the ensemble g(r) bins.
  {
    ngr=0;		    // Number of datasets.
    int number_of_bins = int((bz/2.0)/bin_width);
    cout << "The number of bins is: " << number_of_bins << endl;
    gr.resize(number_of_bins);
    for (int i = 0; i < number_of_bins; i++)
    {
        gr[i].number = i + 1;
        gr[i].width = bin_width;
        gr[i].position = (gr[i].width) * i;           // Assigns bin positions.  See note directly below as well.
        gr[i].population = 0.0;                       // NB moved moving median (center of bins) to post-normalization.
    }
  }
  if (stageFlag == 1)	// Bin all distances contained in the input file.
  {
    // To change between poly- and monodisperse, change lines noted in right margin.
    for (unsigned int i = 0; i < Particle2_List_Bridged.size(); i++)                      // For 1 type, change here (1).
    {
      for (unsigned int j = i+1; j < Particle2_List_Bridged.size(); j++)                    // For 1 type, change here (2).
        {
            VECTOR3D r_vec = Particle2_List_Bridged[i].posvec - Particle2_List_Bridged[j].posvec;     // For 1 type, change here (3).
            if (r_vec.x>bx/2) r_vec.x -= bx;
            if (r_vec.x<-bx/2) r_vec.x += bx;
            if (r_vec.y>by/2) r_vec.y -= by;
            if (r_vec.y<-by/2) r_vec.y += by;
            if (r_vec.z>bz/2) r_vec.z -= bz;
            if (r_vec.z<-bz/2) r_vec.z += bz;
            double r=r_vec.GetMagnitude();
            if (r < bz/2.0)
            {
                int bin_number = ceil((r/bin_width));
                gr[bin_number - 1].population = gr[bin_number - 1].population + 2/BridgedCount;  // For 1 type, optionally change here (3.1).
                //  Added per-addition normalization by "CondensedCount" as the N_condensed  varies with dump step (time).
            }
        }
    }
    cout << "\tDataset binning complete for RDF (select type 2-2)." << endl;
  }
  if (stageFlag == 2)	// Normalize each bin count appropriately (see comments).
  {
    cout << "RDF (select type 2-2) calculation ends, beginning normalization & output." << endl;
    int number_of_bins = int((bz/2.0)/bin_width);
    // Output the RDF:
    ofstream grStream("gr_DD_Bridging_dr=0.0005.out", ios::out);
    for (int b = 0; b < number_of_bins; b++)
    { //  Output moving median abscissae (rMedian) & normalize for g(r) using non-moving median shell volumes binned.
      double r = gr[b].position;
      double vol_bin = (4.0 / 3.0) * pi * ((pow(r + bin_width, 3) - pow(r, 3)));  // Concentric shells' interstitial volume.
      double nid = vol_bin * Particle2_Density;                       // For 1 type, change ParticleJ_Density (4).
      double rMedian = r + 0.5*bin_width;                             // Compute the moving median as final abscissae.
      //gr[b].population = gr[b].population / Particle2_List.size();    // Not necessary, done previously per step.
      gr[b].population = gr[b].population / ngr;                      // Normalize by the number of datasets.
      gr[b].population = gr[b].population / nid;                      // Normalize by the expected number in ideal gas.
      grStream << rMedian << "\t" << gr[b].population << endl;
    }
    grStream.close();
    cout << "\tRDF (select type 2-2) file output complete." << endl;
  }
}

//  Compute MSD m(t): mean square displacement of type 1 particles at a given dumpstep (pushes back to the master MSD):
    //  Unlike in the RDF funcs, output is handled in main; this function simply computes that step, adds to master MSD.
void compute_msd(vector<double> &MSD, vector<PARTICLE> &Particle_List, double bx, double by, double bz)
{
    VECTOR3D COM_vec = VECTOR3D(0,0,0);
    vector<double> totSqDisp_Array; // The list containing all particles' total squared displacements.
    double meanSqDisp = 0.0;

    //  Compute the center of mass at the current dump step:
    for(unsigned int i = 0; i < Particle_List.size(); i++) COM_vec += Particle_List[i].posvec ^ (1.0/Particle_List.size());

    //  Iterate over each particle and push back its displacement at the given dumpstep:
    for(unsigned int i = 0; i < Particle_List.size(); i++)
    {   // Compute the displacement vector, account for periodicity:
        VECTOR3D d_vec = Particle_List[i].posvec - Particle_List[i].initPosVec;// - COM_vec;
        //if (d_vec.x>bx/2) d_vec.x -= bx;
        //if (d_vec.x<-bx/2) d_vec.x += bx;
        //if (d_vec.y>by/2) d_vec.y -= by;
        //if (d_vec.y<-by/2) d_vec.y += by;
        //if (d_vec.z>bz/2) d_vec.z -= bz;
        //if (d_vec.z<-bz/2) d_vec.z += bz;
        double tot_SqDisp = (d_vec.x * d_vec.x) + (d_vec.y * d_vec.y) + (d_vec.z * d_vec.z);

        totSqDisp_Array.push_back(tot_SqDisp);
    }

    // Take the sum to compute the total square displacement, simultaneously taking the mean on the RHS:
    for(unsigned int i = 0; i < totSqDisp_Array.size(); i++)
    {
        meanSqDisp += (totSqDisp_Array[i] / Particle_List.size());
    }

    // Append the number to the final MSD list (y-component of the time series):
    MSD.push_back(meanSqDisp);
}

//  Unused functions:
//  Make movie
void make_movie(int num, vector<PARTICLE>& atom, double bx, double by, double bz, ofstream& outdump)
{
    outdump << "ITEM: TIMESTEP" << endl;
    outdump << num - 1 << endl;
    outdump << "ITEM: NUMBER OF ATOMS" << endl;
    outdump << atom.size() << endl;
    outdump << "ITEM: BOX BOUNDS" << endl;
    outdump << -0.5*bx << "\t" << 0.5*bx << endl;
    outdump << -0.5*by << "\t" << 0.5*by << endl;
    outdump << -0.5*bz << "\t" << 0.5*bz << endl;
    outdump << "ITEM: ATOMS index type x y z" << endl;
    string type;
    for (unsigned int i = 0; i < atom.size(); i++)
    {
        if (atom[i].charge > 0)
            type = "1";
        else
            type = "-1";
        outdump << i+1 << "   " << type << "   " << atom[i].posvec.x << "   " << atom[i].posvec.y << "   " << atom[i].posvec.z << endl;
    }
    return;
}

//  Function to sort by particle type (works with single-digit types only, must be in 2nd column):
bool return_Lesser_Type(string line1, string line2)
{ // Boolean indicator (1) if the particle type of line1 is less than that of line 2:
    return line1[line1.find(" ")+1] < line2[line2.find(" ")+1];
}
