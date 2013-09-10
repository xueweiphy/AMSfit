module InterpolData
use params
use DataLoad
implicit none
   public
   
   private i,j,k,yp1,ypn,temp_y, temp_y2,alpha0Start,res
   private alpha0Step,modeeptemp,pp
   integer i,j,k
   integer modeeptemp
   real(kind=8) yp1,ypn
   real(kind=8) temp_y(stepalphalog) 
   real(kind=8) temp_y2(stepalphalog) 
   real(kind=8):: alpha0Start = 1.9
   real(kind=8):: alpha0Step = 0.1
   real(kind=8) pp,res

   !integer modeep
   !common/dsepdndpiusercom/modeep
   !real*8 alphap,emaxp,eminp,ampdnde,normeng
   !common /dsepdndpiuser4com/alphap,emaxp,eminp,ampdnde, normeng
   !real*8 normprim   ! primary electron normaliztion
   !common/normprimcom/normprim
   !integer nbgmodel
   !real*8 Rsun
   !common/dsnpaxiusercom/Rsun,nbgmodel


   !include 'dscraxicom.h'  ! for rzaxiaddlab, rzaxisuf
   !include 'dsdmdcom.h' ! for nptag
   !include 'dshacom.h'
   !include 'eptest.h'

!!! tabulate primary electron index from 2 to 4
!!! tabulate dm or pulsar spectrum index
   !real(kind=8)::Prim_AMSprotonAlpha(stepalpha,AMSprlnum)
   !real(kind=8)::Prim_AMSprotonAlphay2(stepalpha,AMSprlnum)
   real(kind=8)::Prim_AMSelectronAlpha(stepalpha,AMSelnum)
   real(kind=8)::Prim_AMSelectronAlphay2(stepalpha,AMSelnum)
   real(kind=8)::Prim_AMSpositronAlpha(stepalpha,AMSpolnum)
   real(kind=8)::Prim_AMSpositronAlphay2(stepalpha,AMSpolnum)
   real(kind=8)::Prim_AMSepAlpha(stepalpha,AMSeplnum)
   real(kind=8)::Prim_AMSepAlphay2(stepalpha,AMSeplnum)

   !real(kind=8)::dmps_AMSprotonAlpha(sdimlog,AMSprlnum,stepalphalog)
   !real(kind=8)::dmps_AMSprotonAlphay2(sdimlog,AMSprlnum,stepalphalog)
   real(kind=8)::dmps_AMSelectronAlpha(sdimlog,AMSelnum,stepalphalog)
   real(kind=8)::dmps_AMSelectronAlphay2(sdimlog,AMSelnum,stepalphalog)
   real(kind=8)::dmps_AMSpositronAlpha(sdimlog,AMSpolnum,stepalphalog)
   real(kind=8)::dmps_AMSpositronAlphay2(sdimlog,AMSpolnum,stepalphalog)
   real(kind=8)::dmps_AMSepAlpha(sdimlog,AMSeplnum,stepalphalog)
   real(kind=8)::dmps_AMSepAlphay2(sdimlog,AMSeplnum,stepalphalog)


   !real(kind=8)::Prim_AMSproton(AMSprlnum)
   real(kind=8)::Prim_AMSelectron(AMSelnum)
   real(kind=8)::Prim_AMSpositron(AMSpolnum)
   real(kind=8)::Prim_AMSep(AMSeplnum)


   real(kind=8)::Bge_AMSelectron(AMSelnum)
   real(kind=8)::Bge_AMSpositron(AMSpolnum)
   real(kind=8)::Bge_AMSep(AMSeplnum)
   !real(kind=8)::Sece_AMSelectron(AMSelnum)
   !real(kind=8)::Sece_AMSpositron(AMSpolnum)
   !real(kind=8)::Sece_AMSep(AMSeplnum)

   !real(kind=8):: Sece_AMSelectron(AMSelnum)
   !real(kind=8):: Secpo_AMSelectron(AMSelnum)
   !real(kind=8):: Sece_AMSpositron(AMSpolnum)
   !real(kind=8):: Secpo_AMSpositron(AMSpolnum)
   !real(kind=8):: Sece_AMSep(AMSeplnum)
   !real(kind=8):: Secpo_AMSep(AMSeplnum)


contains

!  Set Primary electron 
! First time Set Primary electron, Please use flagnorm =1 
! to set the primary electron normalization
   !subroutine SetPrimeIni(alphap0, flagnorm)
      !integer flagnorm
      !real(kind=8) alphap0, norm0,pp
      !if ( modeep .eq. 9)
         !epgraxiloadck=.true.
         !epcraxitag='dskra'
         !call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep)
      !else 
         !epgraxiloadck=.false.
      !endif
      !epcraxitag='dskra'
      !modeeptemp = modeep
      !modeep = 9 
      !hamwimp = 1.d4
      !npaxidefault=.false.  !do not link to spherical distribution 
      !nisrf=6
!
      !epspecdef=.false.
      !nptag='ferriere' !tag needed for tabulation of the green function
      !nnptag=8
      !alphap = alphap0
      !if ( flagnorm .eq. 1) then
         !pp = AMSelectronE(21)
         !normprim = AMSelectronF(21)/dsepdndpaxi(Rsun,0.d0,pp,ivopt,how)
      !endif
      !modeep = modeeptemp
   !end subroutine

         

!  AMS proton   
   !subroutine SplinePrimAMSproton()
      !do i = 1, AMSprlnum
         !call spline(prim_alpha,prim_AMSprotonAlpha(:,i),stepalpha,  &
            !0.d0,0.d0, prim_AMSprotonAlphay2(:,i))
      !enddo
   !end subroutine
!
   !subroutine SplineDmpsAMSproton()
      !do j = 1, sdimlog 
         !do i = 1, AMSprlnum
            !temp_y = dmps_AMSprotonAlpha(j,i,:)
            !call spline(alphaloglist,temp_y,stepalphalog,            &
               !0.d0,0.d0,temp_y2)  
            !dmps_AMSprotonAlphay2(j,i,:) = temp_y2
         !enddo
      !enddo
   !end subroutine


   subroutine SplineDmpsAMS()
      call SplineDmpsAMSelectron
      call SplineDmpsAMSpositron
      call SplineDmpsAMSAllelectron
   end subroutine

   subroutine SplintDmpsAMS(alpha0)
      real(kind=8) alpha0(sdimlog)
      call SplintDmpsAMSelectron(alpha0)
      call SplintDmpsAMSpositron(alpha0)
      call SplintDmpsAMSAllelectron(alpha0)
   end subroutine

   subroutine SplinePrimAMS()
      call SplinePrimAMSelectron
      call SplinePrimAMSpositron
      call SplinePrimAMSAllelectron
   end subroutine

   subroutine SplintPrimAMS( alphap0, norm0)
      real(kind=8) alphap0,norm0
      call SplintPrimAMSelectron(alphap0,norm0)
      call SplintPrimAMSpositron(alphap0,norm0)
      call SplintPrimAMSAllelectron(alphap0,norm0)
   end subroutine


! AMS electron
   subroutine SplinePrimAMSelectron()
      do i = 1, AMSelnum
         call spline(prim_alpha,prim_AMSelectronAlpha(:,i),stepalpha,&
            0.d0,0.d0, prim_AMSelectronAlphay2(:,i))
      enddo
   end subroutine

   subroutine SplineDmpsAMSelectron()
      do j = 1, sdimlog 
         do i = 1, AMSelnum
            temp_y = dmps_AMSelectronAlpha(j,i,:)
            call spline(alphaloglist,temp_y,stepalphalog,            &
               0.d0,0.d0,temp_y2)  
            dmps_AMSelectronAlphay2(j,i,:) = temp_y2
         enddo
      enddo
   end subroutine


   subroutine SplintPrimAMSelectron(alphap0,norm0)
      real(kind=8) alphap0,norm0
      do i = 1, AMSelnum
         call splint(prim_alpha,Prim_AMSelectronAlpha(:,i),          &
            Prim_AMSelectronAlphay2(:,i),stepalpha,alphap0,          &
            Prim_AMSelectron(i))
      enddo
      Prim_AMSelectron = norm0 * Prim_AMSelectron
   end subroutine

   subroutine SplintDmpsAMSelectron(alphap0)
      real(kind=8) alphap0(sdimlog)
      AMSelectronSpcV =0.d0
      do j = 1, sdimlog 
         do i = 1, AMSelnum
            temp_y = dmps_AMSelectronAlpha(j,i,:)
            temp_y2 = dmps_AMSelectronAlphay2(j,i,:)
            call splint(alphaloglist,temp_y,temp_y2,stepalphalog,    &
               alphap0(j),res)
            AMSelectronSpcV(j,i) = res
         enddo
      enddo
   end subroutine

! AMS positron 
   subroutine SplinePrimAMSpositron()
      do i = 1, AMSpolnum
         call spline(prim_alpha,prim_AMSpositronAlpha(:,i),stepalpha,&
            0.d0,0.d0, prim_AMSpositronAlphay2(:,i))
      enddo
   end subroutine

   subroutine SplineDmpsAMSpositron()
      do j = 1, sdimlog 
         do i = 1, AMSpolnum
            temp_y = dmps_AMSpositronAlpha(j,i,:)
            call spline(alphaloglist,temp_y,stepalphalog,            &
               0.d0,0.d0,temp_y2)  
            dmps_AMSpositronAlphay2(j,i,:) = temp_y2
         enddo
      enddo
   end subroutine

   subroutine SplintPrimAMSpositron(alphap0,norm0)
      real(kind=8) alphap0,norm0
      do i = 1, AMSpolnum
         call splint(prim_alpha,Prim_AMSpositronAlpha(:,i),          &
            Prim_AMSpositronAlphay2(:,i),stepalpha,alphap0,          &
            Prim_AMSpositron(i))
      enddo
      Prim_AMSpositron = norm0 * Prim_AMSpositron
   end subroutine

   subroutine SplintDmpsAMSpositron(alphap0)
      real(kind=8) alphap0(sdimlog)
      AMSpositronSpcV =0.d0
      do j = 1, sdimlog 
         do i = 1, AMSpolnum
            temp_y = dmps_AMSpositronAlpha(j,i,:)
            temp_y2 = dmps_AMSpositronAlphay2(j,i,:)
            call splint(alphaloglist,temp_y,temp_y2,stepalphalog,    &
               alphap0(j),res)
            AMSpositronSpcV(j,i) = res
         enddo
      enddo
   end subroutine

! AMS electron + positron
   subroutine SplinePrimAMSAllelectron()
      do i = 1, AMSeplnum
         call spline(prim_alpha,prim_AMSepAlpha(:,i),stepalpha,&
            0.d0,0.d0, prim_AMSepAlphay2(:,i))
      enddo
   end subroutine

   subroutine SplineDmpsAMSAllelectron()
      do j = 1, sdimlog 
         do i = 1, AMSeplnum
            temp_y = dmps_AMSepAlpha(j,i,:)
            call spline(alphaloglist,temp_y,stepalphalog,            &
               0.d0,0.d0,temp_y2)  
            dmps_AMSepAlphay2(j,i,:) = temp_y2
         enddo
      enddo
   end subroutine

   subroutine SplintPrimAMSAllelectron(alphap0,norm0)
      real(kind=8) alphap0,norm0
      do i = 1, AMSeplnum
         call splint(prim_alpha,Prim_AMSepAlpha(:,i),                &
            Prim_AMSepAlphay2(:,i),stepalpha,alphap0, Prim_AMSep(i))
      enddo
      Prim_AMSep = norm0 * Prim_AMSep
   end subroutine


   subroutine SplintDmpsAMSAllelectron(alphap0)
      real(kind=8) alphap0(sdimlog)
      AMSepSpcV =0.d0
      do j = 1, sdimlog 
         do i = 1, AMSeplnum
            temp_y  =dmps_AMSepAlpha(j,i,:)
            temp_y2 =dmps_AMSepAlphay2(j,i,:)
            call splint(alphaloglist,temp_y,temp_y2,stepalphalog,    &
               alphap0(j),res)
            AMSepSpcV(j,i) = res
         enddo
      enddo
   end subroutine


end module
