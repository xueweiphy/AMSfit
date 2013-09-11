#include "DataExp.h"

DataExp::DataExp(){

   pathexp="/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/"   ;
}

int DataExp::LoadFermi(string Name){
   string fname;
   fname = pathexp  + Name;
   float** fermidata;
   fermidata= read_file(fname,Frow, 8);

   for( int i = 0; i<Frow; i++)
   {
      Feng[i] = fermidata[2][i];
      //Feng2[i] = Feng[i]*(1.-shift);
      Fflux3[i] = fermidata[3][i];
      Fflux0[i] = Fflux3[i]/pow( Feng[i],3);
      //Fflux32[i] = Fflux3[i]*pow(1.-shift,3);
      Fstat[i] = fermidata[4][i];
      Fsysp[i] = fermidata[5][i];
      Fsysm[i] = fermidata[6][i];
      Ferrorp[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);
      Ferrorp0[i] = Ferrorp[i]/ pow( Feng[i],3);
      Ferrorm[i]=  sqrt( Fsysp[i]*Fsysp[i]+ Fstat[i] * Fstat[i]);
   }

   return 0;
}

int DataExp::LoadAMS(string Name ){
   float** AMSdata;
   string fname;
   fname = pathexp  + Name;
   AMSdata= read_file(fname,AMSrow,5);
   for( int i = 0; i<AMSrow; i++)
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
   return 0;

}




int DataExp::LoadPamelaElectron(string Name ){
   float** pameladata;
   string fname;
   fname = pathexp  + Name;
   pameladata= read_file(fname,pamelarow,7);
   for( int i = 0; i<pamelarow; i++)
   {
      Rmin[i] = pameladata[0][i];
      Rmax[i] = pameladata[1][i];
      pamelaE[i] = pameladata[2][i];
      pamelaflux[i] = pameladata[4][i];
      pamelaerr[i] = sqrt( pow(pameladata[5][i],2)
         + pow(pameladata[6][i],2)  );
   }
   return 0;
}




int DataExp::LoadHess(string Name ){
   float** hessdata;
   string fname;
   fname = pathexp  + Name;
   hessdata= read_file(fname,hessrow,6);
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



   return 0;
}























