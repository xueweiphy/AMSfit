#include "dnde.h"

void dnde() {

//////  Load dnde amplitude data
   string prefix;
   string str_temp;
   prefix = "../chains/dmstep10_";
   str_temp = prefix + "dnde_mod1.dat";
   TNtuple  dndef( "dndef","dNdE","row:amplitude:sigma");
   dndef.ReadFile( str_temp.c_str()); 
    
   int sdim; 
   sdim = dndef.GetEntries();
   int sdimextra;
   int sdimT;
   sdimextra=3;
   sdimT= sdim;
   sdim = sdim-sdimextra;

   cout<<"sdim = " <<sdim<<endl;
   float log10e[sdim];
   float eng[sdim];
   float spectrum[sdim];
   float spectrum2[sdim];
   float  sigma[sdim];
   float errorL[sdim];
   float errorR[sdim];
   float *row ;
   int ii;
   float englogE = 3.6;
   float englogstep = 0.3;
   float etpar[sdimextra];
   float etsig[sdimextra];
   float slp, shift,normpri;
      
   slp = 0.;
   shift =0.;
   normpri =0.;
    
   for ( int i =0; i< sdim; i++)
   {
      ii= sdim -i-1;
      dndef.GetEntry(i);
      row = dndef.GetArgs();
      log10e[i] =  englogE-ii*englogstep;
      eng[i] = pow( 10., log10e[i]);
      cout<<"eng "<<ii<<" "<<eng[i]<<endl;
      spectrum[ii] = row[1];
      cout<<"amp "<<ii<<" "<<spectrum[ii]<<endl;
      sigma[ii] = row[2];
      cout<<"sigma "<<ii<<" "<<sigma[ii]<<endl;
   }

   for ( int i = sdim ; i< sdimT; i++)
   {
      dndef.GetEntry(i);
      row = dndef.GetArgs();
      etpar[i-sdim] = row[1];
      cout<<"addition "<<i<<"  "<< etpar[i-sdim];
      etsig[i-sdim] = row[2];
      cout<<", sigma "<< etsig[i-sdim]<<endl;;
   }
   shift=etpar[0];
   slp=etpar[1];
   normpri = etpar[2];
 

   for ( int i=0; i< sdim; i++)
   {
      spectrum2[i]=0.;
      for ( int j=i; j< sdim; j++)
      {
         spectrum2[i]=spectrum2[i]+spectrum[j];
      }
   }

////// bestfit point   

   float** bestamp;
   str_temp = prefix + "dnde_best.dat";
   bestamp = read_file(str_temp.c_str(),sdimT,2);
   cout<<"===bestfit==="<<endl;
   for ( int i=0; i< sdim; i++)
   {
      ii= sdim -i-1;
      spectrum[ii] = bestamp[1][i];
      cout<<"amp "<<ii<<" "<<spectrum[ii]<<endl;

   }
   shift=bestamp[1][sdim];
   cout<<"shfit = " <<shift<<endl;
   slp=bestamp[1][sdim+1];
   normpri = bestamp[1][sdim+2];


//////  Load Fermi data 
   float** fermidata;
   fermidata= read_file("../eptestF/FERMI_total_12months.dat",42,8);
   //TNtuple  fermi( "dndef","dNdE"
   //,"Emin:Emax:Energy:Flux3:stat:sysp:sysm:smooth");
   //  fermi.Readfile("../eptestF/FERMI_total_12months.dat");
   float Feng[42],Fflux3[42],Fstat[42],Fsysp[42],Fsysm[42],
         Ferrorp[42],Ferrorm[42];
   float Feng2[42], Fflux32[42];
   int Frow;
   Frow=42;
   for( int i = 0; i<42; i++)
   {
      Feng[i] = fermidata[2][i]; 
      Feng2[i] = Feng[i]*(1.-shift);
    
      Fflux3[i] = fermidata[3][i];
      Fflux32[i] = Fflux3[i]*pow(1.-shift,3);
      Fstat[i] = fermidata[4][i];
      Fsysp[i] = fermidata[5][i];
      Fsysm[i] = fermidata[6][i];
      Ferrorp[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);
      Ferrorm[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);
   }

//////  Load AMS data
   float** AMSdata;
   int AMSrow =65;
   AMSdata= read_file("../eptestF/AMS02.dat",65,5);
   float AMSerr1[65],AMSerr2[65],AMSerr[65],AMSpr[65],
         AMSEmin[65],AMSEmax[65],AMSE[65];
   for( int i = 0; i<65; i++)
   {
      AMSEmin[i] = AMSdata[0][i] ;
      AMSEmax[i] = AMSdata[1][i] ;
      AMSE[i] = (AMSEmin[i]+AMSEmax[i])*0.5 ;
      AMSpr[i] = AMSdata[2][i] ;
      AMSerr1[i] = AMSdata[3][i] ;
      AMSerr2[i] = AMSdata[4][i] ;
      AMSerr[i] = sqrt( AMSerr1[i]*AMSerr1[i]+ AMSerr2[i]*AMSerr2[i]);
      //  cout<<AMSerr[i]<<endl;
   }
              
  
//////  Load pamela electron data
   float** pameladata;
   int pamelarow =39;
   pameladata= read_file("../eptestF/PAMELA_electrons_v2.dat",39,7);
   float Rmin[39],Rmax[39],pamelaerr[39],pamelaflux[39],
         pamelaE[39];
   for( int i = 0; i<39; i++)
   {
      Rmin[i] = pameladata[0][i];
      Rmax[i] = pameladata[1][i];
      pamelaE[i] = pameladata[2][i];
      pamelaflux[i] = pameladata[4][i];
      pamelaerr[i] = sqrt( pow(pameladata[5][i],2)
         + pow(pameladata[6][i],2)  );
      // cout<<AMSerr[i]<<endl;
   }
              
  
//////  Load HESS data
   float** hessdata;
   int hessrow =9;
   hessdata= read_file("../eptestF/HESS_v2.dat",9,6);
   float hessE[hessrow],hessflux3[hessrow],
         hesserr1m[hessrow],hesserr1p[hessrow],hesserrm[hessrow],
         hesserr2m[hessrow],hesserr2p[hessrow],hesserrp[hessrow];
             
   for( int i = 0; i<hessrow; i++)
   {
      hessE[i] = hessdata[0][i];
      hessflux3[i] = hessdata[1][i];
      hesserr1m[i] = hessdata[2][i];
      hesserr1p[i] = hessdata[3][i];
      hesserr2p[i] = hessdata[4][i];
      hesserr2m[i] = hessdata[5][i];
      hesserrm[i] = sqrt( pow(hesserr1m[i],2)
         +  pow(hesserr2m[i],2)  );
      hesserrp[i] = sqrt( pow(hesserr1p[i],2)
         +  pow(hesserr2p[i],2)  );
   }
              
  



 
     
//////  Load electron and positron background

   float** kradata;
   int bglnum = 146;
   str_temp = prefix + "spectrum_kra.dat";
   kradata = read_file( str_temp.c_str(),bglnum,3);
   float bgeflux[bglnum],bgpflux[bglnum],bgeng[bglnum];
   float bgepflux3[bglnum];
   for ( int i =0; i< bglnum; i++)
   {
      bgeng[i] = kradata[0][i];
      bgeflux[i] = kradata[1][i]*1.e4;
      bgpflux[i] = kradata[2][i]*1.e4;
      bgepflux3[i] = (bgpflux[i]+bgeflux[i])* pow( bgeng[i],3);
   }

      //TNtuple  bge( "bge","KRA electron","flux");
      //bge.ReadFile( "../eptestF/propmod/electronflux_KRA004_BG.dat"); 
      //TNtuple  bgp( "bgp","KRA positron","flux");
      //bgp.ReadFile( "../eptestF/propmod/positronflux_KRA004_BG.dat"); 
      //TNtuple  bgeng0( "bgeng0","KRA Energy list","eng");
      //bgeng0.ReadFile( "../eptestF/propmod/KRA004_ENG.dat"); 
/*
   for ( int i =0; i< bglnum; i++)
   {
      bgeng0.GetEntry(i);
      row = bgeng0.GetArgs();
      bgeng[i]= row[0];

      bge.GetEntry(i);
      row = bge.GetArgs();
      bgeflux[i]= row[0];
      if( i<10) cout<<setw(4)<<i<< setw(14)<<bgeflux[i]<<endl;
      bgeflux[i]= bgeflux[i]*pow( bgeng[i]/10., slp);
      if( i<10)cout<<setw(4)<<i<< setw(14)<<bgeflux[i]<<setw(8)<<slp<<endl;
          

      bgp.GetEntry(i);
      row = bgp.GetArgs();
      bgpflux[i]= row[0];

      bgepflux3[i] = (bgpflux[i]+bgeflux[i])* pow( bgeng[i],3);
   }
*/     


//////  Load step function flux data
    
   float** stepffile;
   int dmnum;
   dmnum = 146;
   int  stepnum=sdim;
   str_temp = prefix + "spectrum_stepf.dat" ;
   stepffile= read_file(str_temp.c_str(),stepnum,dmnum);
   float DMflux[dmnum];
   float DMflux3[dmnum];
   float DMeng[dmnum];
   float epTotal[dmnum];
   float eTotal[dmnum];
   float epratio[dmnum];
//     float DMengstep; 
//      DMengstep = 0.025;

   for ( int j=0; j<dmnum; j++)
   {  
      DMeng[j] = bgeng[j]; 
      DMflux[j] =0.;
      for ( int i = 0 ; i < sdim; i++)
      {
         ii= sdim -i-1;
         DMflux[j] =  DMflux[j]+spectrum[ii]*stepffile[j][i]*1.e4;
      }
      DMflux3[j] = 2.*DMflux[j]*pow(DMeng[j],3);
      epTotal [j] = DMflux3[j]+ bgepflux3[j];
      eTotal [j] = DMflux[j]+ bgeflux[j];
      epratio [j] =( DMflux[j]+ bgpflux[j]) /
         (2.* DMflux[j]+ bgpflux[j]+  bgeflux[j]) ;
    }

////// Load step function flux at Fermi energy
   float** stepffile2;
   str_temp = prefix + "spectrum_FERMI.dat" ;
   stepffile2= read_file(str_temp.c_str(),
      stepnum,42);
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
  

   TCanvas *mycanvas = new TCanvas("mycanvs","A Simple Graph with error bars",200,10,1200,1800);
    //   mycanvas->SetFillColor(42);
 //  mycanvas->SetGrid();
 //  mycanvas->GetFrame()->SetFillColor(18);
//    mycanvas -> Divide(1,2);
 //   mycanvas -> cd(1);
//  mycanvas->GetFrame()->SetBorderSize(12);
   mycanvas->SetFillColor(18);
TPad*  pad1 = new   TPad("pad1","dnde",     0.02,0.67,0.48,0.99,0);
 TPad*  pad2 = new TPad("pad2","Fermi data",0.02,0.34,0.48,0.66,0);
 TPad*  pad3 = new TPad("pad3","dnde data", 0.52,0.67,0.98,0.99,0);
 TPad*  pad4 = new TPad("pad4","AMS data",  0.52,0.34,0.98,0.66,0);
 TPad*  pad5 =new TPad("pad5","pamela data",0.02,0.01,0.48,0.33,0);
 TPad*  pad6 = new TPad("pad6","Fermi data",0.52,0.01,0.98,0.33,0);
//TH1F *hr = mycanvas->DrawFrame(10,10,4000,100);
//mycanvas->DrawFrame(10,0,4000,200);

   pad1 -> Draw();
   pad2 -> Draw();
   pad3 -> Draw();
   pad4 -> Draw();
   pad5 -> Draw();
   pad6 -> Draw();

//////
////// Pad 1
//////

   pad1->cd();
//   TGraphErrors graph(n_points,x_vals,y_vals,NULL,y_errs);
  
   gROOT -> SetStyle("Plain");
   //mycanvas-> GetPad(1)->SetBorderMode(-1);
   pad1->GetFrame()->SetBorderMode(-1);
   //mycanvas-> GetPad(1)->SetBorderSize(5);
   pad1->GetFrame()->SetBorderSize(5);
   //  mycanvas-> GetPad(1)->SetFillColor(0);
   //  pad1->GetFrame()->SetFillColor(42);
   pad1-> SetLogx();
   pad1-> SetLogy();
   TGraphErrors *gr = new TGraphErrors(sdim,eng,spectrum,NULL,sigma);

   TGraph * dNdEp = new TGraph(sdim,eng, spectrum2);
   gr->SetTitle("dNdE Mode1");
   gr->GetXaxis() -> SetTitle("Energy  [GeV]");
   gr->GetXaxis() -> CenterTitle(1);
   gr->GetYaxis() -> SetTitle("dN/dE  [GeV^{-1}]");
   gr->GetYaxis() -> CenterTitle(1);
//   gr->SetMarkerColor(kBlue);
   gr->SetMarkerColor ( 4 ) ;
   gr->SetMarkerStyle(2);
   gr->SetLineColor ( 4 ) ;
//   gr->SetMarkerStyle(kOpenCircle);
//   gr->SetLineColor ( kBlue ) ;
//   gr->SetLineColor ( 4 ) ;
   gr->SetLineWidth(1);
//   gr->Draw("CLP");
   gr->Draw("ALP");
//    dNdEp -> Draw("C");
  
///////////////////   Pad 2
//////////////////    Fermi + HESS data
  
 //   mycanvas -> cd(2);
   pad2 -> cd();
   TGraphAsymmErrors *fermiP = new TGraphAsymmErrors(Frow,Feng,
          Fflux3,NULL,NULL,Ferrorm,Ferrorp);
   TGraphAsymmErrors *fermiP2 = new TGraphAsymmErrors(Frow,Feng2,
          Fflux32,NULL,NULL,Ferrorm,Ferrorp);
   TGraphAsymmErrors *hessP = new TGraphAsymmErrors(hessrow,hessE,
          hessflux3,NULL,NULL,hesserrm,hesserrp);
   fermiP->SetTitle("Fermi");
   pad2-> SetLogx();
  //pad2-> SetLogy();
   fermiP->GetXaxis() -> SetTitle("Energy  [GeV]");
   fermiP->GetXaxis() -> CenterTitle(1);
   fermiP->GetYaxis() -> SetTitle("#Phi E^{3}  [m^{-2}s^{-1}sr^{-1}GeV^{2}]");
   fermiP2->SetMarkerStyle(2);
   fermiP2->SetLineColor ( 4 ) ;
   fermiP2->SetMarkerColor ( 4 ) ;

   hessP->SetLineColor ( 6 ) ;
   hessP->SetMarkerColor ( 6 ) ;
   hessP->SetMarkerStyle(33);

   fermiP->SetMarkerStyle(2);
   fermiP->SetLineColor ( 2 ) ;
   fermiP->SetMarkerColor ( 2 ) ;
   fermiP->GetYaxis() -> CenterTitle(1);
   TAxis *axis = fermiP -> GetXaxis();
   axis -> SetLimits(10.,5000.);
   fermiP -> GetHistogram()->SetMaximum(250.);
   fermiP -> GetHistogram()->SetMinimum(30.);

//  mycanvas-> GetPad(2)->SetLogy();
    fermiP -> Draw("AP");
   fermiP2 -> Draw("P");
    hessP -> Draw("P");
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
 


///////////////////   Pad 3
//////////////////    step function spectrum

 
    pad3 -> cd();
   float stfluxtemp[dmnum];
            
 TGraph *stepfP0;
 TGraph *stepfP;
  pad3-> SetLogx();
  pad3-> SetLogy();
       for ( int j= 0; j< stepnum; j++){
        for ( int i = 0; i< dmnum; i++) {
            stfluxtemp[i] =stepffile[i][j]*pow(DMeng[i],2);
	    //cout<<stfluxtemp[i]<<"  "<<DMeng[i]<<endl;
            }
   if (j==0)   
         {
 stepfP0= new TGraph(dmnum,DMeng,stfluxtemp );
        stepfP0->Draw("AC");
   stepfP0 -> GetHistogram()->SetMaximum(1.e-4);
   stepfP0 -> GetHistogram()->SetMinimum(1.e-12);
        }
       else
          {
 stepfP= new TGraph(dmnum,DMeng,stfluxtemp );
      stepfP->Draw("C");
          }
        }
  stepfP0->SetTitle("dNdE");
    
    TAxis *axis3 = stepfP0 -> GetXaxis();
    axis3 -> SetLimits(5.,5000.);

   mycanvas->Update();



///////////////////   Pad 4
/////////////////    AMS data
    pad4 -> cd();
  pad4-> SetLogx();
  pad4-> SetLogy();
   TGraphErrors *amsP = new TGraphErrors(AMSrow,AMSE,AMSpr,NULL,AMSerr);

  amsP->SetTitle("AMS");
   amsP->GetXaxis() -> SetTitle("Energy  [GeV]");
   amsP->GetXaxis() -> CenterTitle(1);
   amsP->GetYaxis() -> SetTitle("#Phi (e^{+}) / ( #Phi(e^{-})+ #Phi(e^{+}) )");
   amsP->GetYaxis() -> CenterTitle(1);
   amsP->SetMarkerStyle(2);
   amsP->SetLineColor ( 2 ) ;
   amsP->SetMarkerColor ( 2 ) ;

    amsP -> Draw("AP");
    TAxis *axisams = amsP -> GetXaxis();
    axisams -> SetLimits(5.,5000.);
 TGraph * eprP = new TGraph(bglnum,bgeng, epratio);
    eprP -> Draw("C");

///////////////////   Pad 5
/////////////////    pamela electron data
    pad5 -> cd();
  pad5-> SetLogx();
  pad5-> SetLogy();

   TGraphErrors *pamelaP = new TGraphErrors(pamelarow,pamelaE,
            pamelaflux,NULL,pamelaerr);
  pamelaP->SetTitle("PAMELA electron");
   pamelaP->GetXaxis() -> SetTitle("Energy  [GeV]");
   pamelaP->GetXaxis() -> CenterTitle(1);
   pamelaP->GetYaxis() -> SetTitle("#Phi (e^{-}) ");
   pamelaP->GetYaxis() -> CenterTitle(1);
   pamelaP->SetMarkerStyle(2);
   pamelaP->SetLineColor ( 2 ) ;
   pamelaP->SetMarkerColor ( 2 ) ;

    TAxis *axis5 = pamelaP -> GetXaxis();
    axis5 -> SetLimits(5.,5000.);
//   pamelaP -> GetHistogram()->SetMaximum(250.);
   pamelaP -> GetHistogram()->SetMinimum(1.e-8);
    pamelaP -> Draw("AP");
 TGraph * etotalP = new TGraph(bglnum,bgeng, eTotal);
    etotalP -> Draw("C");



 // Build and Draw a legend 
//      TLegend leg(.1 ,.7 ,.3 ,.9 ,"Lab. Lesson 1"); 
 //    leg.SetFillColor (0) ; 
 //gr->SetFillColor (0) ; 
 //leg.AddEntry(gr,"Exp. Points"); 
 //leg.AddEntry(&f,"Th. Law"); 
 //leg.DrawClone("Same"); 
         // Draw an arrow on the canvas 

   mycanvas->Update();
   mycanvas -> Print("dnde_mod1_fermi.pdf");
}

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

/*
 #ifndef __CINT__ 
 int main(){ 
    dnde() ; 
 } 
 #endif
*/
