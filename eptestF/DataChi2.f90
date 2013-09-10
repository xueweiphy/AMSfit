module DataChi2
use params
use InterpolData
implicit none
public
   private temp1,chi2,i,j,temp2
   real(kind=8):: temp1,chi2,temp2
   integer:: modeep
   integer i,j
   common/dsepdndpiusercom/modeep

contains

   function chi2AMSep()
      real(kind=8):: chi2AMSep
      chi2 = 0.d0
      do i  = AMSeplnum, 1, -1
         temp1 = 0.d0
         if ( AMSepE(i) < EngCut) EXIT
         if ( (modeep .eq.8) .or. ( modeep .eq. 10) )then
            temp1 = temp1 + 2.* AMSepSpcV(1,i)
         else 
            do j  = 1, sdim1
               temp1 =temp1 + 2.* amptemp(j)*AMSepSpcV(j,i)
            enddo
         endif
         temp2 = temp1+Bge_AMSep(i) + Sece_AMSep(i)-AMSepF(i)
         !write(*,*) i, 'flux residual ', temp2, Bge_AMSep(i),AMSepF(i)
         if ( temp2 .gt. AMSepF(i)) then
            chi2 = chi2 + ( temp2/AMSepErru(i))**2
         else
            chi2 = chi2 + ( temp2/AMSepErrd(i))**2
         endif
         !write(*,*) i, 'AMSep chi2 ', chi2
      enddo
      !write(*,*) "chi2AMSep=",chi2
      chi2AMSep =chi2
   end function

end module
