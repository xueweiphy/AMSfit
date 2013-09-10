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


   LoadFermi();
   LoadAMS();
   LoadPamelaElectron();
   LoadHess();

   cout<<"Multinest has totally"<<dimT <<" parameters to be fitted"<<endl;
   cout<<dim<<" of which are DM spectrum amplitudes."<<endl;
   cout<<"the highest energy of the spetrum is 10^"<<englogE;
   cout<<" with logstep "<< Logstep<<endl;



}

int DATA::LoadMod1(string Name){
   string fname;
   fname = path2+prefix + Name; 
   dataMod1 = read_file( fname.c_str(), dimT, 3);
   int ii; 
   for ( int i ; i < sdim;  i++){
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
   
   fname = prefix + Name;
   kradata = read_file( fname.c_str(),bglnum,3);
   for ( int i =0; i< bglnum; i++)
   {
      bgeng[i] = kradata[0][i];
      bgeflux[i] = kradata[1][i]*1.e4;
      bgpflux[i] = kradata[2][i]*1.e4;
      bgepflux3[i] = (bgpflux[i]+bgeflux[i])* pow( bgeng[i],3);
   }
  // Load Step function flux data 
   fname = prefix + Name2; 
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



















   
