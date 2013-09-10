#ifndef _DataExp_
#define _DataExp_
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


using namespace std;

class DataExp{
protected:
   string pathexp;
   static const int Frow = 42;
   double Feng[Frow],Fflux3[Frow],Fstat[Frow],Fsysp[Frow],Fsysm[Frow],
          Ferrorp[Frow],Ferrorm[Frow];
   double Fflux0[Frow],Fflux02[Frow],Ferrorp0[Frow],Ferrorp32[Frow];
   double Feng2[Frow], Fflux32[Frow];
   
   static const int AMSeplnum = 59;
   double AMSepE[AMSeplnum],AMSepF[AMSeplnum],AMSepErru[AMSeplnum],
         AMSepErrd[AMSeplnum];
   
   static const int AMSrow =65;
   float AMSerr1[AMSrow],AMSerr2[AMSrow],AMSerr[AMSrow],AMSpr[AMSrow],
         AMSEmin[AMSrow],AMSEmax[AMSrow],AMSE[AMSrow];
   

   static const int pamelarow =39;
   float Rmin[pamelarow],Rmax[pamelarow],pamelaerr[pamelarow],
         pamelaflux[pamelarow], pamelaE[pamelarow];


   static const int hessrow =9;
   float hessE[hessrow],hessflux3[hessrow],
         hesserr1m[hessrow],hesserr1p[hessrow],hesserrm[hessrow],
         hesserr2m[hessrow],hesserr2p[hessrow],hesserrp[hessrow];




public:
   DataExp();
   int LoadFermi( string Name = "FERMI_total_12months.dat");
   int LoadAMS( string Name = "AMS02.dat");
   int LoadPamelaElectron( string Name = "PAMELA_electrons_v2.dat");
   int LoadHess( string Name = "HESS_v2.dat");
   
};

#endif

