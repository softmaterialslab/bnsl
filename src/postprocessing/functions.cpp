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

/* later dev:
extern vector<PARTICLE> Particle2_List_Condensed, Particle2_List_Bridged;
extern double pi, Na, unitlength, CondensedCount, BridgedCount;
// Quantities used within functions, but that can't be re-declared or initialized each time:
extern vector<double> CondensedCount_VsTime_List, BridgedCount_VsTime_List; // List of the number of selected particles of type 2 vs. used dump step.
*/

extern double pi, Na, unitlength;

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
    vector<string> instDumpLines;             		//  List of all lines to be included in the instantaneous dump file.
    ifstream movieStream;                     		//  Master dumpfile from which to extract instantaneous dump files.
    movieStream.open("outfiles/dump.melt");   		//  Opens the master dumpfile.
    std::string instFileName = "dumpfiles/";
    ofstream instStream(instFileName, ios::out);  	//  The stream that will be output to the all-types file.
    string tempString;                        		//  A temporary string that will be updated as the file is read.
    int j = 0;                                		//  Dummy (line/stream position) variable.
    double fileNumber;                        		//  The file number (if rounded, dump step) being exported.
    int fileNumberInt = 0;

    //  While the movie stream hasn't reached its end, continue importing and exporting requested content:
    while (!movieStream.eof())
    {
       j++;                                                      //  Iterate to track (line) position of the stream.
       fileNumber = ((double)j/(9.0+particleCount) - 1.0);       //  The current dumpstep (net # lines / # per step).
       fileNumberInt = (int) fileNumber;
       getline(movieStream, tempString);
       instDumpLines.push_back(tempString);
       if(((j%(9 + particleCount)) == 0) && (initDumpStep <= fileNumber) && (fileNumber < (initDumpStep + dataSetCount))){
	  instFileName = "dumpfiles/" + std::to_string(fileNumberInt) + ".melt";
       	  instStream.open(instFileName);                                                          // Open stream to file.
       	  sort(instDumpLines.end()-particleCount,instDumpLines.end(),numeric_string_compare);     // Sort by index.
       	  for(vector<string>::iterator it = instDumpLines.begin() ; it != instDumpLines.end(); ++it)
 	     instStream << *it << endl;
       	     instStream.close();
       	  if (fileNumberInt%100 == 0) 
	     cout << "\tParticle-index-sorted time-ensemble sample number " << fileNumberInt <<  endl;
    	}
    	
    	if((j%(9 + particleCount)) == 0) 
	   instDumpLines.clear();
    	if(fileNumber > (initDumpStep + dataSetCount)) 
	   break; // Breaks the loop if the fileNumber exceeds that requested.
     }

     movieStream.close();
     cout << "\nSample generation done. Total number of samples: " << fileNumberInt+1 << "\n" << endl;
}

//  Compute RDF g(r): pair correlation function for same type (self) particles (say, E-E)
void compute_self_gr(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE>& Particle_List, long double bx, long double by, long double bz, double bin_width, long double Particle_Density, ofstream& grStream)
{
  //  Initialize the ensemble g(r) bins
  if (stageFlag == 0){				
    ngr=0;						// Number of datasets
    int number_of_bins = int((bz/2.0)/bin_width);	
    //cout << "The number of bins for g(r) computation are " << number_of_bins << "\n" << endl;
    gr.resize(number_of_bins);
    for (int i = 0; i < number_of_bins; i++)
    {
      gr[i].number = i+1;
      gr[i].width = bin_width;
      gr[i].position = (gr[i].width)*i;		// Assigns bin positions. See note directly below as well
      gr[i].population = 0.0;						// Moving center of bins eval to post-normalization
    }
  }

  //  Bin all distances contained in the input file.
  if (stageFlag == 1){				
  for (unsigned int i = 0; i < Particle_List.size(); i++)
    {
      for (unsigned int j = i+1; j < Particle_List.size(); j++)	
      {
        VECTOR3D r_vec = Particle_List[i].posvec - Particle_List[j].posvec;  
        if (r_vec.x>bx/2) r_vec.x -= bx;
        if (r_vec.x<-bx/2) r_vec.x += bx;
        if (r_vec.y>by/2) r_vec.y -= by;
        if (r_vec.y<-by/2) r_vec.y += by;
        if (r_vec.z>bz/2) r_vec.z -= bz;
        if (r_vec.z<-bz/2) r_vec.z += bz;
        double r=r_vec.GetMagnitude();
        //if (r < bz/2.0 - 1)		// avoiding the g(r) calculation at the largest r
        unsigned int bin_number = ceil((r/bin_width));
	if (bin_number <= gr.size())
           gr[bin_number - 1].population = gr[bin_number - 1].population + 2;  
      }
    }
  }

  //  Normalize each bin count appropriately
  //  Output moving median abscissae (rMedian) & normalize for g(r) using non-moving median shell volumes binned
  if (stageFlag == 2){	
    cout << endl << "RDF (self) calculation ends, beginning normalization & output." << endl;

    int number_of_bins = int((bz/2.0)/bin_width);
    //ofstream grStream("outfiles/gr_VV_dr=0.005.out", ios::out);

    for (int b = 0; b < number_of_bins; b++)
    { 
      double r = gr[b].position;
      double vol_bin = (4.0 / 3.0) * pi * ((pow(r + bin_width, 3) - pow(r, 3)));	// Concentric shells' interstitial volume
      double nid = vol_bin * Particle_Density;                       			// For 1 type, change ParticleJ_Density (4)
      double rMedian = r + 0.5*bin_width;                             			// Compute the moving median as final abscissae
      gr[b].population = gr[b].population / Particle_List.size();    			// Normalize by the number of particles
      gr[b].population = gr[b].population / ngr;                      			// Normalize by the number of datasets
      gr[b].population = gr[b].population / nid;                      			// Normalize by the expected number in ideal gas
      //if (rMedian <= 9)                                             			// Limiting to not show data at artificial cutoff
        grStream << rMedian << "\t" << gr[b].population << endl;
    }
    grStream.close();
    cout << "RDF (self) file output complete." << endl;
  }
}

//  Compute RDF g(r): pair correlation function for diff type (cross) particles (say, virus-dendrimer or E-K)
void compute_cross_gr(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE>& Particle1_List, vector<PARTICLE>& Particle2_List, long double bx, long double by, long double bz, double bin_width, long double Particle1_Density, long double Particle2_Density)
{
  // Initialize the ensemble g(r) bins
  if (stageFlag == 0)	
  {
    ngr=0;		    				// Number of datasets
    int number_of_bins = int((bz/2.0)/bin_width);
    cout << "The number of bins is: " << number_of_bins << endl;
    gr.resize(number_of_bins);
    for (int i = 0; i < number_of_bins; i++)
    {
      gr[i].number = i+1;
      gr[i].width = bin_width;
      gr[i].position = (gr[i].width)*i;             	// Assigns bin positions
      gr[i].population = 0.0;                       
    }
  }
  
  // Bin all distances contained in the input file
  if (stageFlag == 1)	
  { 
    for (unsigned int i = 0; i < Particle1_List.size(); i++)                      
    {
      for (unsigned int j = 0; j < Particle2_List.size(); j++)                    
      {
        VECTOR3D r_vec = Particle1_List[i].posvec - Particle2_List[j].posvec;     
        if (r_vec.x>bx/2) r_vec.x -= bx;
        if (r_vec.x<-bx/2) r_vec.x += bx;
        if (r_vec.y>by/2) r_vec.y -= by;
        if (r_vec.y<-by/2) r_vec.y += by;
        if (r_vec.z>bz/2) r_vec.z -= bz;
        if (r_vec.z<-bz/2) r_vec.z += bz;
        double r=r_vec.GetMagnitude();
        //if (r < bz/2.0 - 1)				// edge leads to some non-trivial errors, avoiding near gr cutoff calc
        unsigned int bin_number = ceil((r/bin_width));
	if (bin_number <= gr.size())
           gr[bin_number - 1].population = gr[bin_number - 1].population + 1;  // all particles in j loop are distinct from ith particle
      }
    }
  }

  // Normalize each bin count appropriately
  if (stageFlag == 2)	
  {
    cout << endl << "RDF (cross) calculation ends, beginning normalization & output." << endl;
    int number_of_bins = int((bz/2.0)/bin_width);
    ofstream grStream("outfiles/gr_EK_dr=0.005.out", ios::out);

    for (int b = 0; b < number_of_bins; b++)
    { 
      double r = gr[b].position;
      double vol_bin = (4.0 / 3.0) * pi * ((pow(r + bin_width, 3) - pow(r, 3)));  	// Concentric shell interstitial volume
      double nid = vol_bin * Particle2_Density;                       			// Relative to type 1 particle, use Particle2_Density
      double rMedian = r + 0.5*bin_width;                             			// Compute the moving median as final abscissae
      gr[b].population = gr[b].population / Particle1_List.size();    			// Normalize by the number of particles (type 1)
      gr[b].population = gr[b].population / ngr;                      			// Normalize by the number of datasets
      gr[b].population = gr[b].population / nid;                      			// Normalize by the expected number in ideal gas
      //if (rMedian <= 9)                                             			// Limiting to not show data at large r
      	 grStream << rMedian << "\t" << gr[b].population << endl;
    }
    grStream.close();
    cout << "RDF (cross) file output complete." << endl << endl;
  }
}

/* not used now; later devel:
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

//  Function to sort by particle type (works with single-digit types only, must be in 2nd column):
bool return_Lesser_Type(string line1, string line2)
{ // Boolean indicator (1) if the particle type of line1 is less than that of line 2:
    return line1[line1.find(" ")+1] < line2[line2.find(" ")+1];
}
*/
