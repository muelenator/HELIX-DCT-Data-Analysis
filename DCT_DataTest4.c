/*
 * DCT_DATATEST4.c
 * 
 * Reads in one event from proto-DCT data and picks out the event on each wire
 * Shows histogram of event (only on wires where it registered)
 *
 */



#include <stdio.h>
#include <stdlib.h>


#define INIT_ROI(X) X = {.minval = 10000, .minloc = -1, .maxval = 10000, .maxloc = -1, .t_eStart = -1, .t_eEnd=0, .spikeOver=0}
#define NUMWIRES 8
#define NUMTSTEPS 1000
#define NUMADCS 32
#define ROISIZE 25


{
	/* Open the data file */
	char		infile[560] = "NI_PDCT_17.txt";
	ifstream	in;
	in.open(infile);
	
	/* Pre-Defines based on data structure */
	int		adc_offsets[NUMADCS] = { -1, 1, -6, -7, 3, 4, -2, -1, 0, 1, 
									 -3, -2, -1, -1, -1, -1 };
	int		threshOffset[NUMWIRES] = {0,-7,2,0,3,2,-1,-7};
	int		thresh[NUMWIRES] = {0};
	int		threshval = -20;
	int 	safeMinimum = -2000;
	int		safeMaximum = 25;
	int 	min_eStart = 2;
	int 	min_eEndVoltage = 4;
	
	/* Saves information per-event. Used for each adc & each wire. 24 used total */
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

	Int_t		adc[NUMADCS][NUMTSTEPS];	// Stores adc readings. Last 16 of each row aren't used. 
	char		cNum[10];					// Used to store file read
	ROI 		ROI_adc[2 * NUMWIRES];		// Stores relevant data of each ADC
	ROI			ROI_sum[NUMWIRES];			// Stores relevant data of each wire, or sum of ADCS
	bool 		waveGood[NUMWIRES];			// Keeps track of events above threshold, but aren't flukes
	
	/* Offset ADC thresholds */
	for (int i=0; i<NUMWIRES; i++) thresh[i] += threshval + threshOffset[i];
	
	/* Temporary: Skip n events. NI_PDCT_17's first registered event is no good */
	int n = 1;
	for (int i = 0 ; i < n*NUMTSTEPS; i++){
		in.ignore(256, '\n');
	}
	
	/* Get data from one event */
	for (int t = 0; t < NUMTSTEPS; t++)
	{
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
				// std::cout << "Wire " << w << " wiresum = " << ROI_sum[w].wireSum[t] << std::endl;
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
		
		/* Readout of where/if an event was registered. t=974 means no event registered */
		std::cout << ROI_sum[w].t_eStart << " " << ROI_sum[w].t_eEnd;
		std::cout << " " << ROI_sum[w].minval << std::endl;

		
	}
	
	/* Plot region of interest */
	TCanvas		*c1 = new TCanvas("c1", "DCT: Canvas 1", 20, 20, 800, 800);
	gStyle->SetOptStat(0);
	c1->Divide(2,4,.01,0.01);

	TH2F *h[NUMWIRES];
	char *histname = new char[10];
	char *titlename = new char[10];
	
	for (int w=0; w<NUMWIRES; w++) {
		if (waveGood[w]) {
			sprintf(histname,"histo%d",w+1);
			sprintf(titlename,"Wire %d",w+1);
			h[w] = new TH2F(histname, titlename, ROISIZE, ROI_sum[w].t_eStart, ROI_sum[w].t_eEnd,
			300, -250, 50);
			
			for (int t=ROI_sum[w].t_eStart; t <= ROI_sum[w].t_eEnd; t++) {
				h[w]->Fill(t, ROI_sum[w].wireSum[t]);
			}
			c1->cd(w+1);
			h[w]->SetDirectory(0);
			h[w]->Draw("BOX");
		}
	}

}