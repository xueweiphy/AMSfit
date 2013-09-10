//∗∗∗∗ Builds a graph with errors, displays it and saves it as image. ∗∗∗ 
 // first , include some header files (within CINT, these will be ignored) 
//#include "TCanvas.h" 
//#include "TROOT.h"
//#include "TGraphErrors.h" 
// #include "TF1.h" 
// #include "TLegend.h" 
// #include "TArrow.h"
// #include "TLatex.h"
 void macro1(){
// The values and the errors on the Y axis
 TCanvas *mycanvas = new TCanvas();
const int  n_points =10;
double x_vals [ n_points]= {1,2,3,4,5,6,7,8,9,10};
double  y_vals [ n_points]= {6 ,12 ,14 ,20 ,22 ,24 ,35 ,45 ,44 ,53};
double y_errs [ n_points]= {5 ,5 ,4.7 ,4.5 ,4.2 ,5.1 ,2.9 ,4.1 ,4.8 ,        5.43};
// Instance of the graph
TGraphErrors *gr = new TGraphErrors ( n_points , x_vals , y_vals , NULL , y_errs ) ; 


  gr->SetTitle("TGraphErrors Example; lengtht  [cm];Arb. Units");
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(kOpenCircle);
  gr -> SetLineColor ( kBlue ) ;
  gr->Draw("ALP");

//graph.SetTitle("Measurement XYZ");
// Make the plot estetically better
//gROOT−>Reset ; 
//gROOT−>SetStyle("Plain") ; 
//graph.SetMarkerStyle ( kOpenCircle ) ; 
//graph.SetMarkerColor ( kBlue ) ;
//graph.SetLineColor ( kBlue ) ;
//	The	canvas	on	which	we ' l l	draw	the	graph
//TCanvas* mycanvas=new TCanvas ();

// Draw the graph !
//graph.DrawClone("APE");

// TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10    ,700,500);
// TCanvas *c1 = new TCanvas();

// Define a linear function
//TF1 f("Linear law" ,"[0]+x*[1]" ,.5 ,10.5) ; 

// Let's make the funcion line nicer 
//f.SetLineColor(kRed);
// f.SetLineStyle(2);
	// Fit it to the graph and draw it 
//	graph.Fit(&f); 
//	f.DrawClone("Same"); 
	// Build and Draw a legend 
//	TLegend leg(.1 ,.7 ,.3 ,.9 ,"Lab. Lesson 1"); 
//	leg.SetFillColor (0) ; 
//	graph.SetFillColor (0) ; 
//	leg.AddEntry(&graph,"Exp. Points"); 
//	leg.AddEntry(&f,"Th. Law"); 
//	leg.DrawClone("Same"); 
	// Draw an arrow on the canvas 
//	TArrow arrow(8,8,6.2,23,0.02,"----|>"); 
//	arrow.SetLineWidth (2) ; 
//	arrow.DrawClone ( ) ; 
	// Add some text to the plot 
	//TLatex text(8.2,7.5,"#splitline{Maximum}{Deviation}"); 
	//text.DrawClone ( ) ; 
	//mycanvas−>Print("graph_with_law.pdf") ; 
 //       c1-> Update();
	mycanvas−>Print("graph_with_law.pdf") ; 
 } 
/*
 #ifndef __CINT__ 
 int main(){ 
    macro1() ; 
 } 
 #endif
*/
