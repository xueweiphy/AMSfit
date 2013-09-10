#include "main.h"

int DrawData(){
   
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
