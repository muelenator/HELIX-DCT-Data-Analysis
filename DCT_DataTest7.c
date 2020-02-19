/*
 * DCT_DATATEST7.c
 *
 * Reads in all events from data files. Finds mins and maxs for each event, both
 * per wire and per event (max/min of all wires). 
 *
 * Plots histogram event start time per wire
 * Plots histogram of drift time per wire
 * Integrates and plots dN/dt (dV/dt) to get time-distance relation
 *
 */

#include <stdio.h>
#include <stdlib.h>

#define INIT_ROI(X)     \
  X = {.minval = 10000, \
       .minloc = -1,    \
       .maxval = 10000, \
       .maxloc = -1,    \
       .t_eStart = -1,  \
       .t_eEnd = 0,     \
       .spikeOver = false}

#define NUMWIRES 8
#define NUMTSTEPS 1000
#define NUMADCS 32
#define ROISIZE 25
#define NUMEVENTS 10000

/*******************************************************************************
 * Saves information per-event. Used for each adc & each wire. 3*NUMWIRES used
 * total. ROI = region of interest, has a max bin size (above)
 * Note: 'minvals' is just the maximum voltages (are packaged negatively).
*******************************************************************************/
typedef struct ROI {
  int minval;              // Wire minimum
  int minloc;              // Wire minimum bin number
  int maxval;              // Wire maximum
  int maxloc;              // Wire maximum bin number
  int t_eStart;            // Event start time
  int t_eEnd;              // Event end time
  bool spikeOver;          // If event ends before ROISIZE, set this to true
  int wireSum[NUMTSTEPS];  // Stores combined voltage of left + right ADCs
} ROI;

/*******************************************************************************
 * Stores information about valid events
*******************************************************************************/
typedef struct per {
  int minvals[NUMEVENTS];   // Event minimum
  int integral[NUMEVENTS];  // Integral of event voltage
  int dn_dt[NUMEVENTS];     // Integral of event voltage time derivative
} per;

/*******************************************************************************
 * Sets several histogram properties
*******************************************************************************/
TH1F * histEditor(int hNum, const char* type, const char* label,
								  const char* axis, int n, int nmin, int nmax) {
	TH1F* h;
  // Naming schemes
  char* histname = new char[10];
  char* titlename = new char[10];
	char* axisname = new char[100];
  sprintf(histname, "%s %d", type, hNum + 1);
  sprintf(titlename, "%s %d", label, hNum + 1);
	sprintf(axisname, "%s", axis);
  // Hist creation + ownership
  h = new TH1F(histname, titlename, n, nmin, nmax);
  h->SetDirectory(0);
  h->GetXaxis()->SetTitle(axisname);
	
	return h;
}

/*******************************************************************************
 * Main
*******************************************************************************/
void DCT_DataTest7(){
  /*****************************************************************************
  * Opens data file
  *****************************************************************************/
  char infile[560] = "NI_PDCT_17.txt";
  ifstream in;
  in.open(infile);

  /*****************************************************************************
  * Pre-defines based on data stucture and PDCT people
  *****************************************************************************/
  int adc_offsets[NUMADCS] = {
      -1, 1, -6, -7, 3,  4,  -2, -1,
      0,  1, -3, -2, -1, -1, -1, -1};  // Offset voltages, came with data set
  int threshOffset[NUMWIRES] = {
      0, -7, 2, 0, 3, 2, -1, -7};  // Threshold offsets, came with data set
  int thresh[NUMWIRES] = {0};      // Stores threshold values for each wires
  int threshval = -50;      // (PARAM) Min voltage to be considered an event
  int safeMinimum = -2000;  // Anything outside the safe min/max gets thrown out
  int safeMaximum = 25;
  int min_eStart = 2;  // # ROI start time = minloc - min_eStart
  int threshFrac = 8;  // Inverse % of threshold for event to be considered over

  /*****************************************************************************
  * Stores information per event
  *****************************************************************************/
  Int_t adc[NUMADCS][NUMTSTEPS];  // Stores adc readings. Last 16 of each row
                                  // aren't used.
  char cNum[10];                  // Used to store file read
  ROI ROI_adc[2 * NUMWIRES];      // Stores relevant data of each ADC
  ROI ROI_sum[NUMWIRES];    // Stores relevant data of each wire, or sum of ADCs
  bool waveGood[NUMWIRES];  // Keeps track of events above threshold, but aren't
                            // flukes

  /*****************************************************************************
  * Stores information about good events
  *****************************************************************************/
  per minPerWire[8];
  per minPerEvent;

  memset(minPerEvent.minvals, 0, sizeof *minPerEvent.minvals);
  memset(minPerEvent.integral, 0, sizeof *minPerEvent.integral);
	memset(minPerEvent.dn_dt, 0, sizeof *minPerEvent.dn_dt);
  for (int w = 0; w < NUMWIRES; w++) {
    memset(minPerWire[w].minvals, 0, sizeof *minPerWire[w].minvals);
    memset(minPerWire[w].integral, 0, sizeof *minPerWire[w].integral);
		memset(minPerWire[w].dn_dt, 0, sizeof *minPerWire[w].dn_dt);
  }

  /*****************************************************************************
  * Sets up histograms
  *****************************************************************************/
  // Canvases
  TCanvas* c1 = new TCanvas("c1", "t_d Start Time Per Wire", 20, 20, 800, 800);
  TCanvas* c2 = new TCanvas("c2", "Drift Time Per Wire", 20, 20, 800, 800);
  TCanvas* c3 = new TCanvas("c3", "Time Distance Relation", 20, 20, 800, 800);
  c1->Divide(2, 4, .01, 0.01);
  c2->Divide(2, 4, .01, 0.01);
  c3->Divide(2, 4, .01, 0.01);
  gStyle->SetOptStat(0);

  // Histograms. Hist number corresponds to canvas
  TH1F* h1[NUMWIRES];
  TH1F* h2[NUMWIRES];
  TH1F* h3[NUMWIRES];

  // Histogram properties
  for (int w = 0; w < NUMWIRES; w++) {
    h1[w] = histEditor(w, "StartTimes", "Wire", "Event Start Time (t)", 
											50, 0, 600);
    h2[w] = histEditor(w, "DriftTimes", "Wire", "Drift Time (t)", 
											30, 0, 60);
    h3[w] = histEditor(w, "Radius", "Wire", "Radius ()", 50, -100, 50);
  }

  /*****************************************************************************
  * Offset ADC thresholds (pre-defines)
  *****************************************************************************/
  for (int i = 0; i < NUMWIRES; i++)
    thresh[i] += threshval + threshOffset[i];

  /*****************************************************************************
  * Starts analysis. Goes through data file one event at a time.
  *****************************************************************************/
  for (int event = 0; event < NUMEVENTS; event++) {
    /* Get one event */
    for (int t = 0; t < NUMTSTEPS; t++) {
      for (int iadc = 0; iadc < NUMADCS - 1; iadc++) {
        in.getline(cNum, 10, ',');
        adc[iadc][t] = atoi(cNum);
        adc[iadc][t] -= adc_offsets[iadc];
      }
      in.getline(cNum, 10);
      adc[31][t] = atoi(cNum);
    }

    /* Find the time of the event + min and max vals */
    for (int w = 0; w < NUMWIRES; w++) {
      int Ladc = 2 * w;      // Left adc reading
      int Radc = 2 * w + 1;  // Right adc reading

      INIT_ROI(ROI_adc[Ladc]);  // Inits. values for algorithm + re-usability
      INIT_ROI(ROI_adc[Radc]);
      INIT_ROI(ROI_sum[w]);

      waveGood[w] = true;  // innocent until proven guilty

      /* Find the min + max voltage for each adc + each wire */
      for (int t = 0; t < NUMTSTEPS; t++) {
        int Lval = adc[Ladc][t];
        int Rval = adc[Radc][t];
        ROI_sum[w].wireSum[t] = Lval + Rval;

        /* Check for malfunction */
        if (Lval < safeMinimum || Rval < safeMinimum) {
          waveGood[w] = false;
          break;
        } else if (Lval > safeMaximum || Rval > safeMaximum) {
          waveGood[w] = false;
          break;
        }
        /* Left ADC */
        if (Lval < ROI_adc[Ladc].minval) {
          ROI_adc[Ladc].minval = Lval;
          ROI_adc[Ladc].minloc = t;
        }
        if (Lval > ROI_adc[Ladc].maxval) {
          ROI_adc[Ladc].maxval = Lval;
          ROI_adc[Ladc].maxloc = t;
        }
        /* Right ADC */
        if (Rval < ROI_adc[Radc].minval) {
          ROI_adc[Radc].minval = Rval;
          ROI_adc[Radc].minloc = t;
        }
        if (Rval > ROI_adc[Radc].maxval) {
          ROI_adc[Radc].maxval = Rval;
          ROI_adc[Radc].maxloc = t;
        }
        /* Together now */
        if (ROI_sum[w].wireSum[t] < ROI_sum[w].minval) {
          ROI_sum[w].minval = ROI_sum[w].wireSum[t];
          ROI_sum[w].minloc = t;
        }
      }

      /* Determines the ROI. If an event is <ROISIZE, sets new stopping point */
      for (int t = 0; t < NUMTSTEPS; t++) {
        if (ROI_sum[w].wireSum[t] < thresh[w] && t <= ROI_sum[w].minloc &&
            ROI_sum[w].t_eStart < 0) {
          if (t < min_eStart)
            ROI_sum[w].t_eStart = 0;
          else
            ROI_sum[w].t_eStart = t - min_eStart;

          ROI_sum[w].t_eEnd = ROI_sum[w].t_eStart + ROISIZE;
        }
        /* Find the bin where the event is pretty much over */
        if (ROI_sum[w].t_eStart > 0 &&
            ROI_sum[w].wireSum[t] > thresh[w]/threshFrac &&
            !ROI_sum[w].spikeOver) {
          ROI_sum[w].spikeOver = true;
          ROI_sum[w].t_eEnd = t;
        }
      }

      /* If no event is found, mark the wave bad */
      if (ROI_sum[w].t_eStart < 0) {
        waveGood[w] = false;
      }
      /* If an event is found, add some data */
      else if (ROI_sum[w].spikeOver && waveGood[w]) {
        for (int t = ROI_sum[w].t_eStart; t < ROI_sum[w].t_eEnd; t++) {
          minPerWire[w].integral[event] += ROI_sum[w].wireSum[t];
          minPerWire[w].dn_dt[event] +=
              ROI_sum[w].wireSum[t] - ROI_sum[w].wireSum[t + 1];
        }
      }

      /* Determine minimum per wire & for all wires in this event */
      if (minPerWire[w].minvals[event] > ROI_sum[w].minval && waveGood[w])
        minPerWire[w].minvals[event] = ROI_sum[w].minval;
      if (minPerEvent.minvals[event] > ROI_sum[w].minval && waveGood[w])
        minPerEvent.minvals[event] = ROI_sum[w].minval;
    }

    /* Add values to histograms */
    for (int w = 0; w < NUMWIRES; w++) {
      if (waveGood[w]) {
        h1[w]->Fill(ROI_sum[w].t_eStart);
        h2[w]->Fill(ROI_sum[w].t_eEnd-ROI_sum[w].t_eStart);
        h3[w]->Fill(minPerWire[w].dn_dt[event]);
      }
    }
  }

  /*****************************************************************************
  * Plot histogram of all minvalues on each wire and of each event
  *****************************************************************************/
  for (int w = 0; w < NUMWIRES; w++) {
    c1->cd(w + 1);
    h1[w]->Draw();
  }
  for (int w = 0; w < NUMWIRES; w++) {
    c2->cd(w + 1);
    h2[w]->Draw();
  }
  for (int w = 0; w < NUMWIRES; w++) {
    c3->cd(w + 1);
    h3[w]->Draw();
  }

  // TFile *hfile = new TFile("DCT_Test5.root","RECREATE","DCT Test 5");

  // Add histogram of minvals of each wire individually

  // Add histogram of minvals of each event

  // TTree *tree = new TTree("RelevantInfo","Relevant information");
  // TBranch *b = tree->Branch ("Total number of events",&numEfound,
  // "numevents/i");
  // TBranch *b1 = tree->Branch ("Events on each wire",&numEperWire,
  // "a/i:b/i:c/i:d/i:e/i:f/i:g/i:h/i");
}
