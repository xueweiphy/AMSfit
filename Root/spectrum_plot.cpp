#include"spectrum_plot.h"

void spc_plot()
{
   gROOT->SetStyle("Plain");
   gStyle -> SetPalette(1);
   TCanvas * canvas = new TCanvas( "canvas1", "canvas",200,10,700,500);
   canvas -> SetLogy() ;
   canvas -> SetLogx() ;
   float** spcdata ;
   float** spcdata_A;
   string prefix = "../chains/dmloglog9_3p3_0p3_";
   int fnum;
   fnum = 42;
   int amsnum;
   amsnum =65;




   spcdata = read_file("spectrum_F_test.dat",4,42);
   spcdata_A = read_file("spectrum_A_test.dat",3,65);


   double Feng[42], Ftest1[42],Ftest2[42],Ftest3[42];
   double AMSE[65], Atest[65],Atest0[65];
   for ( int i = 0; i < 42; i++){
      Feng[i] = spcdata[i][0];
      cout<<Feng[i]<<endl;
      Ftest1[i] = spcdata[i][1];
      Ftest2[i] = spcdata[i][2];
      Ftest3[i] = spcdata[i][3];
   }
   for ( int i = 0; i < 65; i++){
      AMSE[i] = spcdata_A[i][0];
      Atest0[i] = spcdata_A[i][1];
      Atest[i] = spcdata_A[i][2];
   }



   TGraph * spc1 = new TGraph(fnum,Feng, Ftest1);
   TGraph * spc2 = new TGraph(fnum,Feng, Ftest2);
   TGraph * spc3 = new TGraph(fnum,Feng, Ftest3);
   TGraph * spcA0 = new TGraph(amsnum,AMSE, Atest0);
   TGraph * spcA1 = new TGraph(amsnum,AMSE, Atest);
   spc1 -> Draw("AC*");
   //spc2 -> SetLineColor(2);
   //spc2 -> Draw("CP");
   spcA0 -> SetLineColor(2);
   spcA0 -> SetLineWidth(2);
   spcA0 -> Draw( "C");
   //spcA1 -> SetLineColor(2);
   //spcA1 -> SetLineWidth(2);
   //spcA1 -> Draw( "C*");
   
   //spc3 -> Draw("CP.");
   
   canvas -> Print("spc_test.pdf");

}



#ifndef __CINT__
void StandaloneApplication ( int argc , char ** argv ) {
// eventually , evaluate the application parameters argc , argv
// ==>> here the ROOT macroiscalled
spc_plot() ;
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

