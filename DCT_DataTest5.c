/*
 * DCT_DATATEST5.c
 * 
 * Reads in all events from data files. Finds mins and maxs for each event, both
 * per wire and per event (max/min of all wires). Frequency of max voltage 
 * plotted in a histogram in both cases.
 *
 */


#include <stdio.h>
#include <stdlib.h>


#define INIT_ROI(X) X = {.minval = 10000, .minloc = -1, .maxval = 10000, .maxloc = -1, .t_eStart = -1, .t_eEnd=0, .spikeOver=0}

#define NUMWIRES 8
#define NUMTSTEPS 1000
#define NUMADCS 32
#define ROISIZE 25
#define NUMEVENTS 10000


{
	/* Open the data file */
	char		infile[560] = "NI_PDCT_17.txt";
	ifstream	in;
	in.open(infile);
	
	/* Pre-Defines based on data structure */
	int		adc_offsets[NUMADCS] = { -1, 1, -6, -7, 3, 4, -2, -1, 0, 1, -3, -2, 
									-1, -1, -1, -1 }; // Offset voltages, came with data set
	int		threshOffset[NUMWIRES] = {0,-7,2,0,3,2,-1,-7}; // Threshold offsets, came with data set
	int		thresh[NUMWIRES] = {0}; // Stores threshold values for each wires
	int		threshval = -100; // Minimum voltage to be considered an event. Change this however
	int 	safeMinimum = -2000; // Anything outside the safe min/max gets thrown out
	int		safeMaximum = 25;
	int 	min_eStart = 2; // # time bins before max value to start region of interest
	int 	min_eEndVoltage = 4; // Volts above threshold value for the event to be considered over
	
	/* Saves information per-event. Used for each adc & each wire. 3*NUMWIRES used total */
	typedef struct ROI {
		int minval;
		int minloc;
		int maxval;
		int maxloc;
		int t_eStart;
		int t_eEnd;
		int spikeOver;
		int wireSum[NUMTSTEPS];
	} ROI;
	
	/* Stores maximums (minimums in the data set) */
	typedef struct perEvent {
		int minvals[NUMEVENTS];		
	} perEvent;

	Int_t		adc[NUMADCS][NUMTSTEPS];	// Stores adc readings. Last 16 of each row aren't used. 
	char		cNum[10];					// Used to store file read
	ROI 		ROI_adc[2 * NUMWIRES];		// Stores relevant data of each ADC
	ROI			ROI_sum[NUMWIRES];			// Stores relevant data of each wire, or sum of ADCS
	bool 		waveGood[NUMWIRES];			// Keeps track of events above threshold, but aren't flukes
	
	// For minimum histograms
	perEvent minPerWire[8];
	perEvent minPerEvent;
	
	memset(minPerEvent.minvals, 0, sizeof *minPerEvent.minvals);
	for (int w=0; w<NUMWIRES; w++) memset(minPerWire[w].minvals, 0, sizeof *minPerWire[w].minvals);
	
	/* Setup histograms */
	nbins = 50; 	// Number of bins per histogram
	minbin = 0;		// Minimum voltage on hist.
	maxbin = 600;	// Max. voltage on hist
	TCanvas		*c1 = new TCanvas("c1", "DCT: Canvas 1", 20, 20, 800, 800);	// Per-wire canvas
	gStyle->SetOptStat(0);
	c1->Divide(2,4,.01,0.01);
	
	TCanvas		*c2 = new TCanvas("c2", "DCT: Canvas 2", 20, 20, 800, 800); // Per-event canvas
	gStyle->SetOptStat(0);

	TH1F *h[NUMWIRES]; // Individual wire histograms. On the same canvas
	char *histname = new char[10];
	char *titlename = new char[10];
	/* Create & Set basic features of histogram */
	for (int w=0; w<NUMWIRES; w++) {
		sprintf(histname,"histo%d",w+1);
		sprintf(titlename,"Wire %d",w+1);
		h[w] = new TH1F(histname,titlename,nbins,minbin,maxbin);
		h[w]->SetDirectory(0);
		h[w]->GetXaxis()->SetTitle("Max Voltage (V) of on wire (per event)");
	}
	
	TH1F *h1 = new TH1F("EventMins", "Max voltage per event (of all wires)",nbins,minbin,maxbin+400); // Min per event
	h1->SetDirectory(0);
	h1->GetXaxis()->SetTitle("Voltage (V)");
	
	/* Offset ADC thresholds */
	for (int i=0; i<NUMWIRES; i++) thresh[i] += threshval + threshOffset[i];
	
	/* Go through each event in the data file */
	for (int event=0; event<NUMEVENTS; event++) {
		/* Get one event */
		for (int t = 0; t < NUMTSTEPS; t++) {
			for (int iadc = 0; iadc < NUMADCS - 1; iadc++)
			{
				in.getline(cNum, 10, ',');
				adc[iadc][t] = atoi(cNum);
				adc[iadc][t] -= adc_offsets[iadc];
			}
			in.getline(cNum, 10);
			adc[31][t] = atoi(cNum);
		}
	
		/* Find the time of the event + min and max vals */
		for (int w=0; w<NUMWIRES; w++) {
			int Ladc = 2*w;		// Left adc reading
			int Radc = 2*w + 1; // Right adc reading
			
			INIT_ROI(ROI_adc[Ladc]);	// Initializes values for algorithm + re-usability
			INIT_ROI(ROI_adc[Radc]);
			INIT_ROI(ROI_sum[w]);
			
			waveGood[w] = true;

			/* Find the min + max voltage for each adc + each wire */
			for (int t=0; t<NUMTSTEPS; t++) {
				int Lval = adc[Ladc][t];
				int Rval = adc[Radc][t];	
				ROI_sum[w].wireSum[t] = Lval + Rval;
				
				/* Check for malfunction */
				if (Lval < safeMinimum || Rval < safeMinimum) {
					waveGood[w] = false;
					break;
				}
				else if (Lval > safeMaximum || Rval > safeMaximum) {
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

			/* Calls the region of interest (ROI) the bins directly after the minval */
			for (int t=0; t<NUMTSTEPS; t++) {
				if (ROI_sum[w].wireSum[t]<thresh[w] && t<=ROI_sum[w].minloc) {
					if (t < min_eStart) ROI_sum[w].t_eStart = 0;
					else ROI_sum[w].t_eStart = t - min_eStart;
					ROI_sum[w].t_eEnd = ROI_sum[w].t_eStart + ROISIZE;
				}
				/* Find the bin where the event is pretty much over */
				else if (ROI_sum[w].t_eEnd && ROI_sum[w].wireSum[t]>thresh[w]+min_eEndVoltage
						&& !ROI_sum[w].spikeOver) {
					ROI_sum[w].spikeOver = t;
				}
			}
			/* If no event is found, mark the wire */
			if (ROI_sum[w].t_eStart < 0) {
				waveGood[w]=false;
				ROI_sum[w].t_eStart = NUMTSTEPS - ROISIZE - 1;
				ROI_sum[w].t_eEnd = ROI_sum[w].t_eStart + ROISIZE;
			}
			
			if (minPerWire[w].minvals[event]>ROI_sum[w].minval && waveGood[w]) minPerWire[w].minvals[event]=ROI_sum[w].minval;
			if (minPerEvent.minvals[event]>ROI_sum[w].minval && waveGood[w]) minPerEvent.minvals[event]=ROI_sum[w].minval;	
		}
		
		/* Add all minvalues to histograms if the event was good/found on wire w */
		bool eventOK = false;
		for (int w=0; w<NUMWIRES; w++) {
			if (waveGood[w]) {
				h[w]->Fill(-minPerWire[w].minvals[event]);
				eventOK=true;
			}
		}
		if (eventOK) h1->Fill(-minPerEvent.minvals[event]);
			
	}
	/* Plot histogram of all minvalues on each wire and of each event */
	for (int w=0; w<NUMWIRES; w++) {
		c1->cd(w+1);
		h[w]->Draw();
	}
	c2->cd();
	h1->Draw();
	
	
	// TFile *hfile = new TFile("DCT_Test5.root","RECREATE","DCT Test 5");
	
	// Add histogram of minvals of each wire individually
	
	// Add histogram of minvals of each event
	
	// TTree *tree = new TTree("RelevantInfo","Relevant information");
	// TBranch *b = tree->Branch ("Total number of events",&numEfound, "numevents/i");
	// TBranch *b1 = tree->Branch ("Events on each wire",&numEperWire, "a/i:b/i:c/i:d/i:e/i:f/i:g/i:h/i");

}




