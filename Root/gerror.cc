  #include <iostream>
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
#include "TFrame.h"



void gerrors() {
   
   TCanvas *mycanvas = new TCanvas("mycanvs","A Simple Graph with error bars",200,10,700,500);

   mycanvas->SetFillColor(42);
   mycanvas->SetGrid();
   mycanvas->GetFrame()->SetFillColor(21);
  mycanvas->GetFrame()->SetBorderSize(12);

 const int  n_points =10;
 double x_vals [ n_points]= {1,2,3,4,5,6,7,8,9,10};
 double  y_vals [ n_points]= {6 ,12 ,14 ,20 ,22 ,24 ,35 ,45 ,44 ,53};
 double y_errs [ n_points]= {5 ,5 ,4.7 ,4.5 ,4.2 ,5.1,2.9,4.1,4.8,5.43};

   TGraphErrors *gr = new TGraphErrors(n_points,x_vals,y_vals,NULL,y_errs);
//   TGraphErrors graph(n_points,x_vals,y_vals,NULL,y_errs);
  gROOT -> SetStyle("Plain");
   gr->SetTitle("TGraphErrors Example; lengtht  [cm];Arb. Units");
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(kOpenCircle);
   gr->SetLineColor ( kBlue ) ;
   gr->Draw("ALP");
 
 //Define a linear function
TF1 f("Linear law" ,"[0]+x*[1]" ,.5 ,10.5) ; 
 
// Let's make the funcion line nicer 
f.SetLineColor(kRed);
 f.SetLineStyle(2);

 // Fit it to the graph and draw it 
  gr->Fit(&f); 
      f.DrawClone("Same"); 

 // Build and Draw a legend 
      TLegend leg(.1 ,.7 ,.3 ,.9 ,"Lab. Lesson 1"); 
     leg.SetFillColor (0) ; 
    gr->SetFillColor (0) ; 
     leg.AddEntry(gr,"Exp. Points"); 
     leg.AddEntry(&f,"Th. Law"); 
     leg.DrawClone("Same"); 
         // Draw an arrow on the canvas 
       TArrow arrow(8,8,6.2,23,0.02,"----|>"); 
       arrow.SetLineWidth (2) ; 
       arrow.DrawClone ( ) ; 
         // Add some text to the plot 
 TLatex text(8.2,7.5,"#splitline{Maximum}{Deviation}"); 
         text.DrawClone ( ) ; 

   mycanvas->Update();
   mycanvas -> Print("example.pdf");
}
 #ifndef __CINT__ 
 int main(){ 
    gerrors() ; 
 } 
 #endif

