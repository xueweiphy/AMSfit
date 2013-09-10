#ifndef _DATA_
#define _DATA_
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <TSystem.h>
#include <TROOT.h>
#include <vector>
#include <new>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TFrame.h"
#include "TApplication.h"
#include "TNtuple.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "Math/Polynomial.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "TSpectrum.h"
#include "read.h"
#include "bilinear.h"
#include "nr.h"
#include "DataExp.h"


using namespace std;

class DATA: public DataExp{
protected:
   string path,path2,path3;
   string prefix;
   string mod1Name;
   string name; // temporary string
   int dim, dimT, dimEx;
   double *eng;
   double *log10E;
   float **dataMod1,
   int tempInt;
   double *spectrumM;
   double *spectrum;
   double *spectrum2;
   double *sigma;
   double *parEx;
   double *sigEx;
   double shiftM, slpM, normpriM;
   double shift, slp, normpri;
   static const int bglnum =  146;
   float bgeflux[bglnum],bgpflux[bglnum],bgeng[bglnum];
   float bgepflux3[bglnum];
   static const int dmnum =  146;
   float DMflux[dmnum];
   float DMflux3[dmnum];
   float DMeng[dmnum];
   float epTotal[dmnum];
   float eTotal[dmnum];
   float epratio[dmnum];




public:
   DATA(int, int, double, double, string Pref= "dmloglog9_3p3_0p3_");
   int LoadBestFit( string BfName = "dnde_best.dat");
   int GenSpectrum2(int);
   int DealFermi();
   int LoadBackground(string Name = "spectrum_kra.dat", 
         Name2 = "spectrum_stepf.dat")
   
};

#endif

