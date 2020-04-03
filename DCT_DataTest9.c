/*
 * DCT_DATATEST9.c
 *
 * Reads in all events from data files. Finds mins and maxs for each event, both
 * per wire and per event (max/min of all wires). 
 *
 * Plots histogram event start time per wire
 * Plots histogram of drift time per wire
 * Integrates and plots dN/dt (dV/dt) to get time-distance relation
 *
 * Uses wires 3,4,5 (which have similar histograms) to plot new R-t relation
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
void DCT_DataTest9(){
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
  int threshval = -80;      // (PARAM) Min voltage to be considered an event
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
  TCanvas* c1 = new TCanvas("c1", "dN/dt Per Wire with fits", 20, 20, 800, 800);
  TCanvas* c2 = new TCanvas("c2", "r-t Per Wire", 20, 20, 800, 800);
	
	TCanvas* c3 = new TCanvas("c3", "r-t for Wires 3,4,5", 20, 20, 800, 800);
	TCanvas* c4 = new TCanvas("c4", "dN/dt for Wires 3,4,5", 20, 20, 800, 800);
  c1->Divide(2, 4, .01, 0.01);
  c2->Divide(2, 4, .01, 0.01);
  gStyle->SetOptStat(0);

  // Histograms. Hist number corresponds to canvas
  TH1F* h1[NUMWIRES];
  TH1F* h2[NUMWIRES];
	TH1F* h3;
	TH1F* h4;

  // Histogram properties
  for (int w = 0; w < NUMWIRES; w++) {
    h1[w] = histEditor(w, "dN/dt", "Wire", "Drift time (t)", 
                      25, 0, 50);
		h2[w] = histEditor(w, "r-t Relation", "Wire", "Drift Time (t)", 
                      30, 0, 60);
		h2[w]->GetYaxis()->SetTitle("R");					
  }
	
	h3 = histEditor(5, "r-t Relation", "Wires 3-", "Drift Time (t)", 
                      30, 0, 60);
	h4 = histEditor(5, "dN/dt", "Wires 3-", "Drift Time (t)", 
                      25, 0, 50);
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
        h1[w]->Fill(ROI_sum[w].t_eEnd-ROI_sum[w].t_eStart);
      }
    }
		
		// Look for events that occured on the middle 3 wires
		if (waveGood[2] && waveGood[3] && waveGood[4]) {
			for (int i=2; i<5; i++) {
				for (int t=0; t<NUMTSTEPS; t++) {
					h3->Fill(t, h1[i]->Integral(0,t));
					h4->Fill(ROI_sum[i].t_eEnd-ROI_sum[i].t_eStart);
				}
			}
		}
  }
	
  /*****************************************************************************
  * Plot histogram of all minvalues on each wire and of each event
  *****************************************************************************/
	/* Max & min drift times to fit */
	int fitMinVals[NUMWIRES] = {12,9,10,9,9,9,9,14};
	int fitMaxVals[NUMWIRES] = {35,35,33,33,33,33,20,40};
	
	/* Create a bunch of functions and their names to fit to */
	TF1  *gauss[NUMWIRES];
	TF1  *quad[NUMWIRES];
	
	TF1  *cheby[NUMWIRES];
	
	TF1  *gaussD[NUMWIRES];
	TF1  *quadD[NUMWIRES];
	
	char* f1name = new char[10];
	char* f2name = new char[10];
	
	char* f3name = new char[10];
  
	char* f1dname = new char[10];
	char* f2dname = new char[10];
	
	/* Start the fitting process for each wire */
  for (int w = 0; w < NUMWIRES; w++) {
		sprintf(f1name, "Gauss%d", w+1);
		sprintf(f2name, "Pol2%d", w+1);
		
		sprintf(f3name, "Cheby%d", w+1);
		
		sprintf(f1dname, "GaussDer%d", w+1);
		sprintf(f2dname, "Pol2Der%d", w+1);
		
		gauss[w] = new TF1(f1name,"gaus", fitMinVals[w], fitMaxVals[w]);
		quad[w] = new TF1(f2name,"pol 2", fitMinVals[w], fitMaxVals[w]);
		
		gauss[w]->SetLineColor(kRed);
		quad[w]->SetLineColor(kBlue);
		
		/* Draw dN/dt with Gaussian and Quadratic curve fits */
		c1->cd(w + 1);
		h1[w]->Draw();
			
		h1[w]->Fit(gauss[w], "R");
		h1[w]->Fit(quad[w], "R+");
		
		auto derivGauss = [gauss, w] (double *x, double *p) {
			return gauss[w]->Derivative(*x);
		};
		auto derivQuad = [quad, w] (double *x, double *p) {
			return quad[w]->Derivative(*x);
		};
		
		cheby[w] = new TF1(f3name, "cheb5", 3, fitMaxVals[w]);
		
		gaussD[w] = new TF1(f1dname, derivGauss, 3, fitMaxVals[w], 1);
		quadD[w] = new TF1(f2dname, derivQuad, 3, fitMaxVals[w], 1);
		
		cheby[w]->SetLineColor(kBlack);
		
		gaussD[w]->SetLineColor(kRed);
		quadD[w]->SetLineColor(kBlue);
		
		for (int t=0; t<NUMTSTEPS; t++) {
			h2[w]->Fill(t, h1[w]->Integral(0,t));
		}
		
		/* Draw R vs. t with Chebyshev fits (all others don't work b.c. 
		 * we need integrals, not derivatives. Canvas 3 fixes this */
		c2->cd(w+1);
		h2[w]->Draw("P");
		h2[w]->Fit(cheby[w], "R");
		// gaussD[w]->Draw("SAME");
		// quadD[w]->Draw("SAME");
  }
	
	/* Fit dN/dt for wires 3,4,5 */
	TF1 *cheb = new TF1("cheb", "cheb5", 3, fitMaxVals[3]);
	TF1 *gauss1 = new TF1("gaussian", "gaus", fitMinVals[3], fitMaxVals[3]);
	
	c4->cd();
	h4->Draw();
	h4->Fit(gauss1,"R");
	
	/* Fit R vs. t using a chebyshev & a Gaussian integral
	 * Neither looks very good, unfortunately */
	auto iG = [gauss1] (double *x, double *p) {
		return gauss1->Integral(0,*x);
	};
	TF1 *iGauss = new TF1("dgaussian", iG, 3, fitMaxVals[3],1);
	
	c3->cd();
	h3->Draw("P");
	h3->Fit(cheb, "R");
	iGauss->Draw("SAME");

}
