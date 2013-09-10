module like

use params
use utils1
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
      double precision ratio , slp
      double precision fermishiftchi2 ,int_fermishift,a,b
      double precision amschi2, pamelachi2 ,hesschi2
      integer modeep,intflag
      common/dsepdndpiusercom/modeep

        
             
        
         
      TwoPi=6.2831853






       !do i=1,sdim
         !temp(i)=(spriorran(i,2)-spriorran(i,1))*Cube(i)+spriorran(i,1)
         ! change cube(i) from range(0-1) to the range spriorran(1) - spriorran(2)
        !end do





!         slhood = - sdim1 / 2d0 * log( TwoPi )
         slhood = 0.
!           slhood = -1000000000.
          
       !rescaling the parameters in unit hypercube according to the prior    
        do i=1,sdim
         temp(i)=(spriorran(i,2)-spriorran(i,1))*Cube(i)+spriorran(i,1)
        end do
        Cube(1:sdim)=temp(1:sdim)
         amptemp(1:sdim)=Cube(1:sdim)
          slp = Cube(sdim1+2)
        
         if ( modeep .eq. 8) then
            !write(*,*) 'modeep', modeep
             intflag =0
             call spectrumdata(intflag)
         endif
         
        !temp1 =fermishift (0.3)
           
         !do i =1 , 42
         !bg_electron(i) = bg_electron(i)*( FEmid(i)/10.d0)**slp
         !end do

         !do i =1 , 65
          !bg_e_AMS(i) = bg_e_AMS(i)* ( AMSE(i)/10.d0)**slp
         !end do
          
         !do i =1 , 39
          !bg_e_pamela(i) = bg_e_pamela(i)* ( pamelaE(i)/10.d0)**slp
         !end do

         chi2 =0.
         a =-0.4d0
         b= 0.4d0
          chi2 = chi2 + fermishiftchi2(  )
          chi2 = chi2 + amschi2(  )
          chi2 = chi2 + pamelachi2(  )
          chi2 = chi2 + hesschi2(  )

             !Fermi chi2
        !do i = 42,1,-1
         !temp1 =0.  ! temp1 represents the total flux
          !if ( FEmid(i) <20.d0) EXIT          
            !do  j =1, sdim1
             !temp1 = temp1 + Cube(j)* 2.*FspcV(j,i)
             !enddo
        !chi2 = chi2 + (( temp1 +bg_electron(i)+bg_positron(i)- &
         !Fflux(i))/Ferrp(i) )**2
        !enddo
       
        


         !    AMS chi2
        !do i = 65,1,-1
         !temp1 =0.  ! temp1 represents the total flux
          !if ( AMSE(i) <20.d0) EXIT          
            !do  j =1, sdim1
             !temp1 = temp1 + Cube(j)* AMSspcV(j,i)
             !enddo
            !ratio = ( temp1 + bg_p_AMS(i)) / ( 2.*temp1 +    &
         !bg_e_AMS(i)*(AMSE(i)/10.d0)**slp + bg_p_AMS(i)  ) 
        !chi2 = chi2 +( (ratio- AMSpr(i))/ AMSerr(i) )**2
        !enddo


         
         !    pamela electron chi2
        !do i = 39,1,-1
         !temp1 =0.  ! temp1 represents the total flux
          !if ( pamelaE(i) <20.d0) EXIT          
            !do  j =1, sdim1
             !temp1 = temp1 + Cube(j)* pamelaspcV(j,i)
             !enddo
        !chi2 = chi2 + (( temp1 +    &
          !bg_e_pamela(i)* ( pamelaE(i)/10.d0)**slp- &
         !pamelaflux(i))/pamelaerr(i) )**2
        !enddo





          
         slhood = slhood - chi2/2.d0
             

end subroutine slikelihood

      
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








           chi2Function=chi2

  end function 


!double precision function fermishift( delta)
 
       !implicit none
       !double precision delta
   ! 



!end function




!=======================================================================

end module like
