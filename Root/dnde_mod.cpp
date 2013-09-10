#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <TSystem.h> 
#include <TROOT.h>
#include <vector>
#include <math.h>
#include "TCanvas.h" 
#include "TGraphErrors.h" 
#include "TF1.h" 
#include "TH1.h" 
#include "TFile.h" 
#include "TLegend.h" 
#include "TArrow.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFrame.h"
#include "TApplication.h"
#include "TNtuple.h"
#include "TMath.h"
#include "read.cpp"
 using namespace std;

void dnde() {


////////////////  Load Fermi data 
      float** fermidata;
 fermidata= read_file("../eptestF/FERMI_total_12months.dat",42,8);
   //TNtuple  fermi( "dndef","dNdE"
     //,"Emin:Emax:Energy:Flux3:stat:sysp:sysm:smooth");
    
//      fermi.Readfile("../eptestF/FERMI_total_12months.dat");
    float Feng[42],Fflux3[42],Fstat[42],Fsysp[42],Fsysm[42],
          Ferrorp[42],Ferrorm[42];
     float *row ;
      int Frow;
        Frow=42;
  for( int i = 0; i<42; i++)
   {
   Feng[i] = fermidata[2][i]; 
   Fflux3[i] = fermidata[3][i];
   Fstat[i] = fermidata[4][i];
   Fsysp[i] = fermidata[5][i];
   Fsysm[i] = fermidata[6][i];
   Ferrorp[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);
   Ferrorm[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);

   }

////////////////  Load dnde amplitude data

   TNtuple  dndef( "dndef","dNdE","row:amplitude:sigma");
    dndef.ReadFile( "../chains/dnde_mod3.dat"); 
    
    int sdim; 
     sdim = dndef.GetEntries();
       cout<<"sdim = " <<sdim<<endl;
   float log10e[sdim];
    float eng[sdim];
    float spectrum[sdim];
     float  sigma[sdim];
    int ii;
    
    for ( int i =0; i< sdim; i++)
   {
     ii= sdim -i-1;
    dndef.GetEntry(i);
    row = dndef.GetArgs();
    log10e[i] =  3.6-ii*0.2;
     
     eng[i] = pow( 10., log10e[i]);
//  cout<<"eng "<<eng[i]<<endl;
     spectrum[ii] = row[1];
  cout<<"amp "<<ii<<" "<<spectrum[ii]<<endl;
     sigma[ii] = row[2];
  cout<<"sigma "<<ii<<" "<<sigma[ii]<<endl;
   }
    
     
////////////////  Load electron and positron background
   TNtuple  bge( "bge","KRA electron","flux");
    bge.ReadFile( "../eptestF/propmod/electronflux_KRA004_BG.dat"); 
   TNtuple  bgp( "bgp","KRA positron","flux");
    bgp.ReadFile( "../eptestF/propmod/positronflux_KRA004_BG.dat"); 
   TNtuple  bgeng0( "bgeng0","KRA Energy list","eng");
    bgeng0.ReadFile( "../eptestF/propmod/KRA004_ENG.dat"); 
     int bglnum = 146;
      float bgeflux[bglnum],bgpflux[bglnum],bgeng[bglnum];
      float bgepflux3[bglnum];
    for ( int i =0; i< bglnum; i++)
     {
       bge.GetEntry(i);
       row = bge.GetArgs();
       bgeflux[i]= row[0];

       bgp.GetEntry(i);
       row = bgp.GetArgs();
       bgpflux[i]= row[0];

       bgeng0.GetEntry(i);
       row = bgeng0.GetArgs();
       bgeng[i]= row[0];
     bgepflux3[i] = (bgpflux[i]+bgeflux[i])* pow( bgeng[i],3);
     }
     


////////////////  Load step function flux data
    
      float** stepffile;
         int dmnum;
            dmnum = 123;
 stepffile= read_file("../chains/spectrum_stepf.dat",20,dmnum);
     float DMflux[dmnum];
     float DMflux3[dmnum];
     float DMeng[dmnum];
      float epTotal[dmnum];
//     float DMengstep; 
//        DMengstep = 0.025;
       for ( int j=0; j<dmnum; j++)
    {  
 //     DMeng[j] = pow( 10., DMengstep*(j-1));
      DMeng[j] = bgeng[j]; 
//         cout<< DMeng[j]<<" "<<endl;
       DMflux[j] =0.;
      for ( int i = 0 ; i < sdim; i++)
         {
        ii= sdim -i-1;
       DMflux[j] =  DMflux[j]+spectrum[ii]*stepffile[j][i]*1.e4;
         }
        DMflux3[j] = 2.*DMflux[j]*pow(DMeng[j],3);
        epTotal [j] = DMflux3[j]+ bgepflux3[j];
    }

///////////////// Load step function flux at Fermi energy
      float** stepffile2;
 stepffile2= read_file("../chains/spectrum_FERMI.dat",20,42);
     float DMflux3_fermi[42];
       for ( int j=0; j<42; j++)
       {
         
       DMflux3_fermi[j] =0.;
      for ( int i = 0 ; i < sdim; i++)
         {
        ii= sdim -i-1;
       DMflux3_fermi[j] =  DMflux3_fermi[j]+
       2.*spectrum[ii]*stepffile2[j][i]*1.e4*pow(Feng[j],3) ;
         }
         }
  

   TCanvas *mycanvas = new TCanvas("mycanvs","A Simple Graph with error bars",200,10,600,1200);
    //   mycanvas->SetFillColor(42);
 //  mycanvas->SetGrid();
 //  mycanvas->GetFrame()->SetFillColor(18);
//    mycanvas -> Divide(1,2);
 //   mycanvas -> cd(1);
//  mycanvas->GetFrame()->SetBorderSize(12);
   mycanvas->SetFillColor(18);
TPad*  pad1 = new TPad("pad1","dnde Fermi",0.05,0.50,0.95,0.95,0);
TPad*  pad2 = new TPad("pad2","Fermi data",0.05,0.05,0.95,0.45,0);
//TH1F *hr = mycanvas->DrawFrame(10,10,4000,100);
//mycanvas->DrawFrame(10,0,4000,200);

   pad1 -> Draw();
   pad2 -> Draw();

///////////////////   Pad 1
//////////////////

    pad1->cd();
   TGraphErrors *gr = new TGraphErrors(sdim,eng,spectrum,NULL,sigma);
//   TGraphErrors graph(n_points,x_vals,y_vals,NULL,y_errs);
  
  gROOT -> SetStyle("Plain");
   gr->SetTitle("Fermi Mode1");
  //mycanvas-> GetPad(1)->SetBorderMode(-1);
     pad1->GetFrame()->SetBorderMode(-1);
  //mycanvas-> GetPad(1)->SetBorderSize(5);
     pad1->GetFrame()->SetBorderSize(5);
//  mycanvas-> GetPad(1)->SetFillColor(0);
//     pad1->GetFrame()->SetFillColor(42);
  pad1-> SetLogx();
  pad1-> SetLogy();
   gr->GetXaxis() -> SetTitle("Energy  [GeV]");
   gr->GetXaxis() -> CenterTitle(1);
   gr->GetYaxis() -> SetTitle("dN/dE  [GeV^{-1}]");
   gr->GetYaxis() -> CenterTitle(1);
//   gr->SetMarkerColor(kBlue);
   gr->SetMarkerColor ( 2 ) ;
   gr->SetMarkerStyle(21);
   gr->SetLineColor ( 2 ) ;
   gr->SetMarkerStyle(kOpenCircle);
//   gr->SetLineColor ( kBlue ) ;
   gr->SetLineColor ( 4 ) ;
   gr->SetLineWidth(1);
//   gr->Draw("CLP");
   gr->Draw("ALP");
  
///////////////////   Pad 2
//////////////////    Fermi Data
  
 //   mycanvas -> cd(2);
    pad2 -> cd();
   TGraphAsymmErrors *fermiP = new TGraphAsymmErrors(Frow,Feng,
          Fflux3,NULL,NULL,Ferrorm,Ferrorp);
  fermiP->SetTitle("Fermi");
  pad2-> SetLogx();
   fermiP->GetXaxis() -> SetTitle("Energy  [GeV]");
   fermiP->GetXaxis() -> CenterTitle(1);
   fermiP->GetYaxis() -> SetTitle("#Phi E^{3}  [m^{-2}s^{-1}sr^{-1}GeV^{2}]");
   fermiP->SetMarkerStyle(2);
   fermiP->SetLineColor ( 2 ) ;
   fermiP->SetMarkerColor ( 2 ) ;
   fermiP->GetYaxis() -> CenterTitle(1);
  //  TAxis *axis = fermiP -> GetXaxis();
  //  axis -> SetLimits(10.,1000.);
   fermiP -> GetHistogram()->SetMaximum(250.);
   fermiP -> GetHistogram()->SetMinimum(30.);

//  mycanvas-> GetPad(2)->SetLogy();
    fermiP -> Draw("AP");
 TGraph * bgepP = new TGraph(bglnum,bgeng, bgepflux3);
    bgepP -> Draw("C");
 TGraph * DMepP = new TGraph(dmnum,DMeng, DMflux3);
    DMepP -> Draw("C");
 TGraph * TepP = new TGraph(dmnum,DMeng, epTotal);
   TepP->SetLineColor ( 2 ) ;
    TepP -> Draw("C");

    //TGraph * DMepP2 = new TGraph(42,Feng, DMflux3_fermi);
    //DMepP2 -> Draw("C");

   mycanvas->Update();
 
 


 // Build and Draw a legend 
//      TLegend leg(.1 ,.7 ,.3 ,.9 ,"Lab. Lesson 1"); 
 //    leg.SetFillColor (0) ; 
 //gr->SetFillColor (0) ; 
 //leg.AddEntry(gr,"Exp. Points"); 
 //leg.AddEntry(&f,"Th. Law"); 
 //leg.DrawClone("Same"); 
         // Draw an arrow on the canvas 

   mycanvas->Update();
   mycanvas -> Print("dnde_mod3_fermi.pdf");
}
/*
#ifndef __CINT__
void StandaloneApplication ( int argc , char ** argv ) {
// eventually , evaluate the application parameters argc , argv
// ==>> here the ROOT macroiscalled
dnde() ;
}

// This is the standard "main" of C++ starting a ROOT application
int main ( int argc , char ** argv ) {
gROOT-> Reset() ;
 TApplication app ( "Root Application" , &argc , argv ) ;
StandaloneApplication ( app.Argc(), app.Argv() ) ;
app.Run () ;
exit(0);
 return 0 ;
 }
 #endif

*/
 #ifndef __CINT__ 
 int main(){ 
    dnde() ; 
 } 
 #endif
