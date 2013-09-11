#include "main.h"

int DrawData(){
   DATA * run1 = new DATA( 8, 11, 3.0, 0.3, "DmWoFermiLoglog8_3p0_0p3_");
   run1 -> DataIni();
   run1 -> RootIni(2,2);
   run1 -> PlotElectronPositron(1);
   run1 -> PlotElectron(2);
   run1 -> PlotRatio(3);
   run1 -> PlotdNdE(4);
   run1 -> PrintCanvas("DmWoFermiLoglog8_3p0_0p3.pdf");
   
   return 0;
}

#ifndef __CINT__
void StandaloneApplication ( int argc , char ** argv ) {
// eventually , evaluate the application parameters argc , argv
// ==>> here the ROOT macroiscalled
   DrawData();

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

/*
 #ifndef __CINT__ 
 int main(){ 
    dnde() ; 
 } 
 #endif
*/
