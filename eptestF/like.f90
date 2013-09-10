module like

use params
use utils1
use DataChi2
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
      implicit none
      
      double precision Cube(nest_nPar),slhood
      double precision temp(sdim),dist,loclik
      integer i,j
      double precision TwoPi
      double precision chi2 
      double precision temp1
      double precision ratio , slp,normpri0
      double precision fermishiftchi2 ,int_fermishift,a,b
      double precision amschi2, pamelachi2 ,hesschi2
      double precision alphalog0(sdimlog)
      integer modeep,intflag,normflag,lockflag
      integer flagoutrange
      common/dsepdndpiusercom/modeep
             
      TwoPi=6.2831853
      normflag =1

!     slhood = - sdim1 / 2d0 * log( TwoPi )
      slhood = 0.
          
!rescaling the parameters in unit hypercube according to the prior    

      do i=1,sdim
         temp(i)=(spriorran(i,2)-spriorran(i,1))*Cube(i)+spriorran(i,1)
      end do

      Cube(1:sdim)=temp(1:sdim)
      amptemp(1:sdim)=Cube(1:sdim)
      slp = Cube(sdim1+2)      ! prim electron injection
      normpri0 = Cube(sdim1+3) ! and its normalization 

              
      if ( (modeep .eq. 3) .or. ( modeep .eq. 4)   &
        .or. (modeep .eq. 10) .or. ( modeep .eq. 11)&
        .or. (modeep .eq. 12)) then
!!         call setprie_exp( slp,normpri0,normflag)

         call splint_prim(slp,normpri0)

!!!!!!!!!!! COMMENT OUT the following line, if you are run AMS new data
!!!         ( NOT include ratio? ???)
!         call SplintPrimAMS(slp,normpri0)
         call setsecep_exp
         call addprisece
         intflag =0
         if ( modeep .eq. 10) then
            call spectrumdm_data(intflag,1)
         endif
         if ( (modeep .eq. 11) .or. (modeep .eq. 12)   ) then
            call alpha_loglog(amptemp(1:sdim1),alphalog0,flagoutrange)
            if ( flagoutrange .eq. 1) then
               slhood  = -1.d100
               return
            else 
            call splint_spcloglog(amptemp(1:sdim1),alphalog0)
            call SplintDmpsAMS(alphalog0)
            endif
                  
         endif

      else if ( modeep .eq. 8) then
            !write(*,*) 'modeep', modeep
         intflag =0
         call spectrumdata(intflag)
      endif

      chi2 =0.
      a =-0.4d0
      b= 0.4d0
      !chi2 = chi2 + fermishiftchi2(  ) 
      chi2 = chi2 + amschi2(  )        
      chi2 = chi2 + pamelachi2(  )    
      !chi2 = chi2 + hesschi2(  )       
      !chi2 = chi2 + chi2AMSep()
   



          
      slhood = slhood - chi2/2.d0
             

end subroutine slikelihood



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function chi2Function(amp)
 implicit none
   double precision amp(nest_nPar)
   double precision  temp1,chi2,ratio
   double precision fermishiftchi2,amschi2,pamelachi2
   double precision int_fermishift ,hesschi2
      
   integer i,j
        
           !call spline(FEmid,Ferrp,flnum,0.d0,0.d0,Ferry2)
          !write(*,*) ' int fermshift FUNCTION',               &
              !int_fermishift( -0.4d0,0.4d0)
          
            
   chi2 =0.
   chi2 = chi2 + fermishiftchi2()
   write(*,*) 'fermi chi2 = ', chi2
   chi2 = chi2 + amschi2()
   write(*,*) 'ams chi2 = ', amschi2()
   chi2 = chi2 + pamelachi2()
   write(*,*) 'pamela chi2 = ', pamelachi2()
   chi2 = chi2 + hesschi2()
   write(*,*) 'hess chi2 = ', hesschi2()
   chi2 = chi2 + chi2AMSep()
   write(*,*) 'AMSep chi2 = ', hesschi2()








           chi2Function=chi2

  end function 


!double precision function fermishift( delta)
 
       !implicit none
       !double precision delta
    



!end function




!=======================================================================

end module like
