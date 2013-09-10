#include "DATA.h"

DATA::DATA(int Sdim, int SdimT ,double LogE, 
         double Logstep, string Pref) :DataExp()
{
   path="/scratch/xuewei/study/Dark/MultiNest_v2.18/"   ;
   path2="/scratch/xuewei/study/Dark/MultiNest_v2.18/chains/"   ;
   path3="/scratch/xuewei/study/Dark/MultiNest_v2.18/Root/"   ;
   prefix = Pref;
   mod1Name = "dnde_mod1.dat";
   dim = Sdim;
   dimT = SdimT ;
   dimEx = dimT - dim;
   englogE = LogE;
   englogStep = Logstep;
   eng = new double [logstep];
   log10E = new double [logstep];
   sigma = new double [dim];
   spectrum = new double [dim];
   spectrumM = new double [dim];
   parEx = new double [dimEx];
   sigEx = new double [dimEx];
   spectrum2 = new double [dim];
   
   int ii;
   for ( int i = 0; i < dim; i++){
      ii = sdim  - i -1;
      log10E[i] =  englogE - ii * logstep
      eng[i] = pow( 10., log10E);
      cout<< "Energy "<< setw(5)<<ii<< setw(10)<<eng[i]<<endl;
   }


   LoadFermi();
   LoadAMS();
   LoadPamelaElectron();
   LoadHess();

   cout<<"Multinest has totally"<<dimT <<" parameters to be fitted"<<endl;
   cout<<dim<<" of which are DM spectrum amplitudes."<<endl;
   cout<<"the highest energy of the spetrum is 10^"<<englogE;
   cout<<" with logstep "<< logstep<<endl;



}

int DATA::LoadMod1(){
   name = path2+prefix + mod1Name; 
   dataMod1 = read_file( name.c_str(), &tempInt, 3);
   int ii; 
   for ( int i ; i < sdim; i ; i++){
      ii = sdim - i -1 ;
      spectrumM[ii] = dataMod1[1][i];
      cout<<"amp"<<setw(5)<<ii<<setw(12)<<spectrum[ii];
      sigma[ii] = dataMod1[2][i];
      cout<<"  ,sigma"<<setw(12)<<sigma[ii]<<endl;
   }

   for ( int i = sdim ; i< sdimT; i++)
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
   name = path2+prefix + BfName;
   bestamp = read_file(str_temp.c_str(),sdimT,2);
   cout<<"===bestfit==="<<endl;
   for ( int i=0; i< sdim; i++)
   {
      ii= sdim -i-1;
      spectrumB[ii] = bestamp[1][i];
      cout<<"amp"<<setw(5)<<ii<<setw(12)<<spectrum[ii]<<endl;

   }
   shiftB=bestamp[1][sdim];
   cout<<"shfit = " <<shiftB<<endl;
   slpB=bestamp[1][sdim+1];
   cout<<"slope = " <<slp<<endl;
   normpriB = bestamp[1][sdim+2];
   cout<<"primary electron normarlization = " <<normpri<<endl;
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

int DATA::LoadEpBackground(string Name){





}















   
