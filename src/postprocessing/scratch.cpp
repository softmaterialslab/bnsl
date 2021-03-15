//  Compute RDF g(r): pair correlation function for select type 2-2 particles (say, condensed-only dendrimer-dendrimer):
void compute_gr_select_22(int stageFlag, vector<BINCONTAINER>& gr, unsigned int ngr, vector<PARTICLE> &Particle1_List, vector<PARTICLE> &Particle2_List, long double bx, long double by, long double bz, double bin_width, long double Particle1_Density, long double Particle2_Density, int initDumpStep)
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

/*
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
*/