#include "DATA.h"

DATA::DATA(int Sdim, int SdimT ,double LogE, 
         double Logstep, string Pref) :DataExp()
{
   path="/scratch/xuewei/study/Dark/MultiNest_v2.18/"   ;
   path2="/scratch/xuewei/study/Dark/MultiNest_v2.18/chains/"   ;
   path3="/scratch/xuewei/study/Dark/MultiNest_v2.18/Root/"   ;
   prefix = Pref;
   //mod1Name = "dnde_mod1.dat";
   dim = Sdim;
   sdim =dim;
   dimT = SdimT ;
   dimEx = dimT - dim;
   double englogE,englogStep;
   englogE = LogE;
   
   englogStep = Logstep;
   eng = new double [dim];
   log10E = new double [dim];
   sigma = new double [dim];
   spectrum = new double [dim];
   spectrumM = new double [dim];
   parEx = new double [dimEx];
   sigEx = new double [dimEx];
   spectrum2 = new double [dim];
   
   int ii;
   for ( int i = 0; i < dim; i++){
      ii = dim  - i -1;
      log10E[i] =  englogE - ii * Logstep;
      eng[i] = pow( 10, log10E[i]);
      cout<< "Energy "<< setw(5)<<ii<< setw(10)<<eng[i]<<endl;
   }



   cout<<"Multinest has totally"<<dimT <<" parameters to be fitted"<<endl;
   cout<<dim<<" of which are DM spectrum amplitudes."<<endl;
   cout<<"the highest energy of the spetrum is 10^"<<englogE;
   cout<<" with logstep "<< Logstep<<endl;



}

int DATA::DataIni(){
   LoadFermi();
   LoadAMS();
   LoadPamelaElectron();
   LoadHess();
   LoadMod1();
   LoadBestFit();
   DealFermi();
   LoadBackground();
   return 0;
}

int DATA::RootIni( int a , int b){
// set canvas number : a * b
   gROOT-> SetStyle("Plain");
   gStyle->SetPadTickY(1);
   gStyle->SetPadTickX(1);
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   gStyle -> SetTitleOffset(0.95,"x");
   gStyle -> SetTitleOffset(0.95,"y");
   gStyle->SetStatStyle(0);
   gStyle->SetTitleStyle(0);
   gStyle->SetTitleBorderSize(0);
   gROOT->ForceStyle();
   gStyle -> SetLegendBorderSize(0);
   gStyle->SetTitleSize(0.05,"x");
   gStyle->SetTitleSize(0.05,"y");
   gStyle -> SetTitleX( 0.1f);
   gStyle -> SetTitleW( 0.8f);
   double ww,wh;
   ww  = 400*a;
   wh  = 300*b;
   padnumx = a;
   padnumy = b;
   //c_cont = new TCanvas("c_cont"," ",1,1);
   c1 = new TCanvas("c1"," ",ww,wh);
   c1 -> Divide(a,b);
   return 0;
}


int DATA::SetPad( TPad* padG){
   padG->SetLeftMargin(0.10);
   padG->SetRightMargin(0.02);
   padG->SetBottomMargin(0.10);
   padG->SetTopMargin(0.08);
   padG -> SetLogx();
   padG -> SetLogy();
   return 0;
}




int DATA::LoadMod1(string Name){
   string fname;
   fname = path2+prefix + Name; 
   dataMod1 = read_file( fname.c_str(), dimT, 3);
   int ii; 
      cout<<"===== Load Mod1 ===="<<endl;
   for ( int i=0 ; i < sdim;  i++){
      ii = sdim - i -1 ;
      spectrumM[ii] = dataMod1[1][i];
      cout<<"amp"<<setw(5)<<ii<<setw(12)<<spectrum[ii];
      sigma[ii] = dataMod1[2][i];
      cout<<"  ,sigma"<<setw(12)<<sigma[ii]<<endl;
   }

   for ( int i = sdim ; i< dimT; i++)
   {
      parEx[i-sdim] = dataMod1[1][i];
      sigEx[i-sdim] = dataMod1[2][i];
      cout<<"extra"<<setw(5)<<i-sdim<<setw(12)<<parEx[i-sdim];
      cout<<" ,sigma"<<setw(12)<<sigEx[i-sdim]<<endl;
   }
   shiftM=parEx[0];
   slpM=parEx[1];
   normpriM = parEx[2];
   return 0;
}

int DATA::GenSpectrum2(int BM){
   // BM :  =1 Best  =0 Mid
   for ( int i=0; i< sdim; i++)
   {
      spectrum2[i]=0.;
      for ( int j=i; j< sdim; j++)
      {
         if ( BM ==0)
            spectrum2[i]=spectrum2[i]+spectrumM[j];
         else
            spectrum2[i]=spectrum2[i]+spectrum[j];
      }
   }

   return 0;
}

int DATA::LoadBestFit(string BfName){
   float** bestamp;
   string fname;
   fname = path2+prefix + BfName;
   bestamp = read_file(fname.c_str(),dimT,2);
   cout<<"===bestfit==="<<endl;
   int ii;
   for ( int i=0; i< sdim; i++)
   {
      ii= sdim -i-1;
      spectrum[ii] = bestamp[1][i];
      cout<<"amp"<<setw(5)<<ii<<setw(12)<<spectrum[ii]<<endl;

   }
   shift=bestamp[1][sdim];
   cout<<"shfit = " <<shift<<endl;
   slp=bestamp[1][sdim+1];
   cout<<"slope = " <<slp<<endl;
   normpri = bestamp[1][sdim+2];
   cout<<"primary electron normarlization = " <<normpri<<endl;
   return 0;
}


int DATA::DealFermi(){
   for( int i = 0; i<Frow; i++)
   {
      Feng2[i] = Feng[i]*(1.-shift);
      Fflux32[i] = Fflux3[i]*pow(1.-shift,3);
   }
   ROOT::Math::Interpolator * interF =
          new ROOT::Math::Interpolator(Frow);
   ROOT::Math::Interpolator *interErr =
          new ROOT::Math::Interpolator(Frow);
   interF -> SetData( Frow, Feng, Fflux0);
   interErr -> SetData( Frow, Feng, Ferrorp0);
   for (int  i = 0; i< Frow; i++){
      Feng2[i] =  Feng[i]*(1.+shift);
      Fflux02[i] = ( interF-> Eval( Feng2[i]));
      Ferrorp32[i] = ( interErr-> Eval( Feng2[i]));
      Fflux32[i] = Fflux02[i]* pow( Feng[i],3);
      Ferrorp32[i] = Ferrorp32[i]* pow( Feng[i],3);
      //cout<<" Ferrorp32  "<<   Ferrorp32[i] <<endl;
      Feng2[i] = Feng[i];
   }
   return 0;
}

int DATA::LoadBackground(string Name, string Name2){
   float** kradata;
   string fname;
   
   fname = path2+ prefix + Name;
   kradata = read_file( fname.c_str(),bglnum,3);
   for ( int i =0; i< bglnum; i++)
   {
      bgeng[i] = kradata[0][i];
      bgeflux[i] = kradata[1][i]*1.e4;
      bgpflux[i] = kradata[2][i]*1.e4;
      bgepflux3[i] = (bgpflux[i]+bgeflux[i])* pow( bgeng[i],3);
   }
  // Load Step function flux data 
   fname =path2+ prefix + Name2; 
   float** stepffile;
   stepffile= read_file(fname.c_str(),sdim,dmnum);
   int ii;
   for ( int j=0; j<dmnum; j++)
   {
      DMeng[j] = bgeng[j];
      DMflux[j] =0.;
      
      for ( int i = 0 ; i < sdim; i++)
      {
         ii= sdim -i-1;
         DMflux[j] =  DMflux[j]+spectrum[ii]*stepffile[j][i]*1.e4;
         //cout<<setw(5)<<i << setw(15)<<spectrum[ii]<<endl;
      }
      DMflux3[j] = 2.*DMflux[j]*pow(DMeng[j],3);
      epTotal [j] = DMflux3[j]+ bgepflux3[j];
      eTotal [j] = DMflux[j]+ bgeflux[j];
      epratio [j] =( DMflux[j]+ bgpflux[j]) /
         (2.* DMflux[j]+ bgpflux[j]+  bgeflux[j]) ;
    }

   return 0; 
}



int DATA::PlotElectronPositron(int pnum, int expnum,string Title){
   padData = ( TPad *) c1 -> GetPad( pnum);
   SetPad(padData);
   padData -> SetLogy(0);
   //padData -> SetLogx();
   c1 -> cd (pnum);
   TGraphAsymmErrors *fermiP = new TGraphAsymmErrors(Frow,Feng,
          Fflux3,NULL,NULL,Ferrorm,Ferrorp);
   //TGraphAsymmErrors *fermiP2 = new TGraphAsymmErrors(Frow,Feng2,
   //Fflux32,NULL,NULL,Ferrorp32,Ferrorp32);
   //TGraphAsymmErrors *hessP = new TGraphAsymmErrors(hessrow,hessE,
   //hessflux3,NULL,NULL,hesserrm,hesserrp);
   //TGraphAsymmErrors *AMSepP = new TGraphAsymmErrors(AMSeplnum,AMSepE,
   //AMSepF,NULL,NULL,AMSepErrd,AMSepErru);
   fermiP->SetTitle(Title.c_str());
   //fermiP-> CenterTitle(1);
   fermiP->GetXaxis() -> SetTitle("Energy  [ GeV ]");
   fermiP->GetXaxis() -> CenterTitle(1);
   fermiP->GetYaxis() -> SetTitle("#Phi E^{3}  [ GeV^{2} m^{-2}s^{-1}sr^{-1} ]");
   fermiP->GetYaxis() -> CenterTitle(1);
   TAxis *axis = fermiP -> GetXaxis();
   axis -> SetLimits(6.,3000.);
   fermiP -> GetHistogram()->SetMaximum(250.);
   fermiP -> GetHistogram()->SetMinimum(30.);
   fermiP -> SetMarkerColor(2);
   fermiP -> SetLineColor(2);
   fermiP -> SetMarkerStyle(26);
   fermiP -> SetMarkerSize(0.5);
   fermiP -> Draw("AP");

   TGraph * DMepP = new TGraph(dmnum,DMeng, DMflux3);
   DMepP -> Draw("C");
   DMepP -> SetLineWidth(1.5);
   
   TGraph * TepP = new TGraph(dmnum,DMeng, epTotal);
   TepP->SetLineColor ( 4 ) ;
   TepP -> SetLineWidth(1.8);
   TepP -> Draw("C");


   return 0;
   
}

int DATA::PlotElectron(int pnum, int expnum,string Title){
   padData = ( TPad *) c1 -> GetPad( pnum);
   SetPad(padData);
   //padData -> SetLogy(0);
   //padData -> SetLogx(0);
   c1 -> cd (pnum);
   TGraphErrors *pamelaP = new TGraphErrors(pamelarow,pamelaE,
            pamelaflux,NULL,pamelaerr);
   pamelaP->SetTitle(Title.c_str());
   pamelaP->GetXaxis() -> SetTitle("Energy  [ GeV ]");
   pamelaP->GetXaxis() -> CenterTitle(1);
   //pamelaP->GetYaxis() -> SetTitle("#Phi (e^{-}) ");
   pamelaP->GetYaxis() -> SetTitle("#Phi [ GeV^{-1} m^{-2}s^{-1}sr^{-1} ]");
   pamelaP->GetYaxis() -> CenterTitle(1);
   TAxis *axis = pamelaP -> GetXaxis();
   axis -> SetLimits(5.,5000.);
   //   pamelaP -> GetHistogram()->SetMaximum(250.);
   pamelaP -> GetHistogram()->SetMinimum(1.e-8);
   pamelaP -> SetMarkerColor(2);
   pamelaP -> SetLineColor(2);
   pamelaP -> SetMarkerStyle(26);
   pamelaP -> SetMarkerSize(0.5);
   pamelaP -> Draw("AP");


   TGraph * etotalP = new TGraph(bglnum,bgeng, eTotal);
   etotalP -> SetLineWidth(1.8);
   etotalP -> Draw("C");

   return 0;
   
}


int DATA::PlotRatio(int pnum, int expnum,string Title){
   padData = ( TPad *) c1 -> GetPad( pnum);
   SetPad(padData);
   padData -> SetLogy(0);
   //padData -> SetLogx(0);
   c1 -> cd (pnum);

   TGraphErrors *amsP = new TGraphErrors(AMSrow,AMSE,AMSpr,NULL,AMSerr);

   amsP->SetTitle(Title.c_str());
   amsP->GetXaxis() -> SetTitle("Energy  [ GeV ]");
   amsP->GetXaxis() -> CenterTitle(1);
   amsP->GetYaxis() -> SetTitle("#Phi (e^{+}) / [ #Phi(e^{-})+ #Phi(e^{+}) ]"    );
   amsP->GetYaxis() -> CenterTitle(1);
   amsP -> SetMarkerColor(2);
   amsP -> SetLineColor(2);
   amsP -> SetMarkerStyle(24);
   amsP -> SetMarkerSize(0.5);
   amsP -> Draw("AP");
   TAxis *axis = amsP -> GetXaxis();
   axis -> SetLimits(5.,1000.);
   //amsP -> GetHistogram()->SetMaximum(250.);
   //amsP -> GetHistogram()->SetMinimum(30.);

   TGraph * eprP = new TGraph(bglnum,bgeng, epratio);
   eprP -> Draw("C");

   return 0;
   
}


int DATA::PlotdNdE(int pnum, string Title){
   padData = ( TPad *) c1 -> GetPad( pnum);
   SetPad(padData);
   //padData -> SetLogy(0);
   //padData -> SetLogx(0);
   c1 -> cd (pnum);
   //for ( int i = 0; i< sdim ; i++){
   //cout<< setw(10)<< eng[i];
   //cout<< setw(10)<< spectrumM[i];
   //cout<< setw(10)<< sigma[i];
   //cout<< endl;
   //}
   TGraphErrors *gr = new TGraphErrors(sdim,eng,spectrumM,NULL,sigma);
   gr->SetTitle(Title.c_str());
   gr->GetXaxis() -> SetTitle("Energy  [ GeV ]");
   gr->GetXaxis() -> CenterTitle(1);
   gr->GetYaxis() -> SetTitle("dN/dE  [ GeV^{-1} ]");
   gr->GetYaxis() -> CenterTitle(1);
   gr->SetMarkerColor ( 4 ) ;
   gr->SetMarkerStyle(7);
   gr->SetLineColor ( 4 ) ;
   gr->SetLineWidth ( 1.6 ) ;
   gr->Draw("AP");
   
   return 0;
}





int DATA::PrintCanvas(const char* name){
   c1 -> Print(name);
   return 0;
}


