



      include 'user-dummy-eptest.h'

      include 'user-not-dummy-eptest.h'

      program main
        use nestwrapper
        use params

ccc routine to find the N-points fit to Fermi, HESS, AMS, PAMELA
ccc
      implicit none
      character*10 locallabel
      integer i,j
ccc
      include 'dscraxicom.h'  ! for rzaxiaddlab, rzaxisuf
      include 'dsdmdcom.h' ! for nptag
      include 'dshacom.h'

c      real*8 dsepdphidphaaxi ! defined in eptest.h
       include 'eptest.h'
      real*8 ppin  ! momentum

      real*8 radius,res1,res2,dsrhosph_user

      real*8 dummyres

      real*8 mchiin,sigmavin,ppup,pplow,resmat(4,100),dsepdndpi
      integer k,chin,chinvec(4),nbins
      real*8 mstepvec(4),resvec(4,100),engstep
      real*8 massvec(31),spcvec(31,146),engvec(121)
    

ccc
ccc for the user-defined ep initial spectrum
ccc
      integer modeep
      common/dsepdndpiusercom/modeep

ccc
      integer nbinsep
      real*8 deltalogS(100),deltalogp(100),logSendpoint,logpendpoint
      common/dsepdndpiuser2com/deltalogS,deltalogp,logSendpoint,
     &  logpendpoint,nbinsep
ccc  for step funtion spectrum
      real*8 mstep1
      common /massstep/ mstep1

ccc   fermi data
       integer flnum
       parameter ( flnum = 42)
         !integer stepnum
         !parameter(stepnum=40)
         !real*8 Famp(stepnum)

       !parameter ( stepnum = 40)
      real*8 FEmin(42), FEmax(42), Fflux0(42),
     & Fstat(42),Fsysp(42), Fsysm(42), Fsmooth(42) ,
     &  fluxtemp,feps, prec

      
c        FEmid(42),Ferrp(42), Ferrm(42),FspcV(20,42),
c     &  Fflux(42),


cccccccccccc   for  AMS data
          integer amsnum
          parameter ( amsnum= 65)
          real*8 AMSEmin(amsnum), AMSEmax(amsnum),
     &    AMSerr1(amsnum),AMSerr2(amsnum)
          real*8 ratiob
c          real*8  bg_p_AMS(amsnum), bg_e_AMS(amsnum)

c         integer flagdata1, flagdata2 
           integer Ferr
!ccc   KRA backgroud data
       integer klnum 
       parameter ( klnum = 146)
       real*8 kra_eng(klnum),kra_positron(klnum),
     & kra_electron(klnum)
c     defined in params.f90  
c      real*8 bg_electron(flnum), bg_positron(flnum)
        real*8 y2e(klnum)
        real*8 y2p(klnum)
         integer m
         real*8 dsepflux,dsf_int

ccc      step function amplitude
          real*8 stepamp(sdim1),stepampmax(sdim1)
          real*8  amptemp1, amptemp2
       real*8 ampbest(sdim), res_sigma(sdim),cube(sdim)
        integer number1
        character(len=3) tempchar
c         real*8 chi2Function

cccccccccc     for pamela data 
          integer pamelanum
          parameter ( pamelanum= 39)
        real*8 Rmin(39),Rmax(39),pamelaEve(39),
     &  pamelaflux0(39),pamelaerr1(39),
     &  pamelaerr2(39)  

cccccccccccc    for HESS data
           integer hessnum
          parameter ( hessnum= 9)
        real*8 hessflux3(9),hesserr1m(9),hesserr1p(9),
     &         hesserr2m(9),hesserr2p(9)

!c initialization of the dark SUSY code
      call dsinit

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc    Read PAMELA data
      open (unit = 114, file = 
     & 'eptestF/PAMELA_electrons.dat',
     & form ='formatted')
       read(114,*)
       do i = 1, 39
       read(114,*) Rmin(i),Rmax(i),pamelaE(i),pamelaEve(i),
     &  pamelaflux0(i), pamelaerr1(i), pamelaerr2(i)
         pamelaflux(i)= pamelaflux0(i)/1.d4
         pamelaerr(i) = dsqrt( pamelaerr1(i)**2 + pamelaerr2(i)**2)
       enddo
       close(114)  


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  Read Fermi data
!c     open ( unit= 4, file = "FERMI_total_12months.dat")
      open (unit = 115, file = 
cc     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/FERMI_total_12months.dat',
     &   form='formatted')
c      read (115,*)
      do i = 1, 42
       read( 115,*)   FEmin(i), FEmax(i), FEmid(i), Fflux0(i),
     & Fstat(i),Fsysp(i), Fsysm(i), Fsmooth(i) 
        Fflux(i)  = Fflux0(i)/ ( FEmid(i)**3) /1.d4 
       Ferrp(i)=dsqrt(Fstat(i)**2+Fsysp(i)**2)/(FEmid(i)**3)/1.d4
       Ferrm(i)=dsqrt(Fstat(i)**2+Fsysm(i)**2)/(FEmid(i)**3)/1.d4
          enddo
       close(115)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  Read AMS data
      open (unit = 119, file = 
cc     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/AMS02.dat',
     &   form='formatted')
      do i = 1, amsnum
       read( 119,*)   AMSEmin(i), AMSEmax(i), AMSpr(i),
     & AMSerr1(i),AMSerr2(i) 
       AMSE(i) = 0.5* ( AMSEmin(i) + AMSEmax(i));
       AMSerr(i)=dsqrt(AMSerr1(i)**2+AMSerr2(i)**2);
          enddo
       close(119)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  Read HESS data
      open (unit = 124, file = 
cc     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/HESS_official.dat',
     &   form='formatted')
          read( 124,*)
      do i = 1, hessnum
       read( 124,*)   hessE(i), hessflux3(i), hesserr1m(i),
     & hesserr1p(i),hesserr2p(i), hesserr2m(i)
        hessflux(i) = hessflux3(i)/(hessE(i)**3)/1.d4
       
       hesserr(i)=dsqrt(hesserr1p(i)**2+hesserr2p(i)**2);
          enddo
       close(124)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  Read Background Model
      open (unit = 125, file = 
c     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/propmod/KRA004_ENG.dat' , form='formatted')
      open (unit = 126, file = 
c     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/propmod/electronflux_KRA004_BG.dat',   form='formatted')
      open (unit = 127, file = 
c     & '/scratch/xuewei/study/Dark/MultiNest_v2.18/eptestF/
     &'eptestF/propmod/positronflux_KRA004_BG.dat',   form='formatted')
        do i =1, klnum
         read(125,*) kra_eng(i)
         write(*,*) i,kra_eng(i)
         read(126,*) kra_electron(i)
         read(127,*) kra_positron(i)
               kra_electron(i) =  kra_electron(i) / 1.d4 
               kra_positron(i) =  kra_positron(i) / 1.d4
!c         write(*,*) kra_eng(i),kra_electron(i)
        enddo
       close(125)
       close(126)
       close(127)
ccc     background flux at Fermi energy 
         write(*,*) 'background electron and positron flux'
         write(*,*) 'energy, bg_electron, bg_positron, Fermi flux'
         call spline(kra_eng,kra_electron,klnum,0.d0,0.d0,y2e)
         call spline(kra_eng,kra_positron,klnum,0.d0,0.d0,y2p)
         do i = 1, flnum
         call splint(kra_eng,kra_electron,y2e,klnum,
     &  FEmid(i),bg_electron(i))
         
         call splint(kra_eng,kra_positron,y2p,klnum,
     &  FEmid(i),bg_positron(i))
            
  
c         write(*,*) FEmid(i),bg_electron(i), bg_positron(i), Fflux(i),
c     &    Fflux(i) - bg_electron(i)-bg_positron(i)
         enddo

cccccccccc   background flux at AMS energy
         do i = 1, amsnum
         call splint(kra_eng,kra_electron,y2e,klnum,
     &  AMSE(i),bg_e_AMS(i))
         
         call splint(kra_eng,kra_positron,y2p,klnum,
     &  AMSE(i),bg_p_AMS(i))
         write(*,*) bg_e_AMS(i), bg_p_AMS(i)
           enddo 
  

cccccccccc   background flux at PAMELA energy
         do i = 1, pamelanum
         call splint(kra_eng,kra_electron,y2e,klnum,
     &  pamelaE(i),bg_e_pamela(i))
           enddo 

cccccccccc   background flux at hess energy
         do i = 1, hessnum
         call splint(kra_eng,kra_electron,y2e,klnum,
     &  hessE(i),bg_e_hess(i))
         call splint(kra_eng,kra_positron,y2p,klnum,
     &  hessE(i),bg_p_hess(i))
           enddo 



  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc halo model:
ccc
      npaxidefault=.true.  ! link to spherical distribution
      npsphdefault=.true.  ! link to rho^2
      rhodefault=.false.   ! do not link to the default distribution but to
                           ! but to dsnpsph_user

!ccc
!ccc choose a halo profile and make a dummy call to it set the 
!ccc galactocentric distance and the inner truncation radius which are
!ccc hardcoded in dsrhosph_user but can be used as global variables
!ccc
      nptag='pbnfw'   ! this is the reference NFW profile used in the
                      ! pbar paper, see user-not-dummy-ep/dsrhosph_user.f
      nnptag=5
      dummyres=dsrhosph_user(1.d0)

!ccc
!ccc choose the diffusion coefficient and halo height 
!ccc

      include 'propmod/pbkra.f'  ! reference kraichan model in the pbar paper
      diffrcep=0.001d0    ! inner radial cut in ep routines in kpc, this value
                          ! should not matter

      epcraxitag='dskra'
      call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep) ! add label for 
                          ! book keeping and correct link to tabulated green
                          ! function

!ccc
!ccc model "H" for energy loss parameters in the mean field model
!ccc
      ivopt=1  ! the matching between  pp and the variable v is done 
               ! assuming fully general momentum loss rate and spatial
               ! diffusion coefficient, with a numerical integral + a
               ! tabulation involved
!ccc for ivopt=1:
!ccc
!ccc some mean from moskalenko talk or porter paper
!ccc
      starlightmean=1.d0-0.25d0 ! optical + IR, eV cm^-3
!ccc
!ccc some 'local' mean in 8*dexp(-(r-r0)/50kpc) *dexp(-|z|/3kpc)
!ccc random + regular, strong talk 
!ccc
      bfieldmean=6.d0      ! \muG

!ccc
!ccc diffusion time for a few values of the pbar kinetic energy:
!ccc
      how=4       ! green function read from table if it exists
      phiin=0.5d0   ! forced field solar modulation parameter 

!ccc
!ccc
!ccc step function for the input spectrum:
!ccc
      epspecdef=.false.
ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccc       change modeep 
ccccccccc   3--- step function
cccccccc    4--- theta function
      modeep=4

      hasv=3.d-26 ! cm^3 s^-1 -> you still need this because dsepdphidphaaxi
                  ! is taken to be proportional to hasv/hamwimp**2
      phiin=0.d0   ! look at the interstellar

!c       habr(15)=0.d0  ! this sets the monochromatic term at E=hamwimp,
                     ! for habr(15)=1.d0 the spectrum is just bypassed
       


ccc
cccc tabulate the step size function. mass starting from few TeV to 
cccc  the mass smaller 10 GeV
ccc   energy from  1 GeV to 1 TeV
ccc
       open (unit = 2, file = 'chains/spectrum_stepf.dat',
     &  form='formatted' , status='Replace')



        write(*,*) ' write spectrum to spectrum_stepf'
    
c      msteptb =0.30d0  ! log mass step
      engstep=0.025d0  ! log energy step



      mstep1 = msteptb
      do k =1, stepnum 
      massvec(k) = 10.d0**( hlog10- msteptb*(k-1))

      hamwimp = massvec(k) 
          
       
c         if ( k==1)   engvec(i) =10.d0**(  engstep*(i-1))
         engstep = 0.01
      do i =1, 146
c       if ( k==1)   engvec(i) =10.d0**(  engstep*(i-1))
c        ppin= engvec(i)
        ppin= kra_eng(i)
       
c        if ( k.eq.1)  write(*,*) 'momentum', ppin,kra_eng(i)

          if ( ppin > massvec(k) ) then
             spcvec(k,i) =0.
            else 
       spcvec(k,i) = dsepdphidphaaxi(ppin,phiin,how)
           end if
       enddo 
      enddo
      DO k = 1 ,stepnum
       WRITE(2,22) (spcvec(k,i),i=1,146)  ! output to file
       enddo
22     FORMAT(146(1x,e13.8))
      close(2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc tabulate the step size function. mass from 1 TeV to 10 GeV
ccc   energy are at Fermi data points
ccc
 
       feps =1.d-3
       prec = 1.d-2
       
       open (unit = 4, file = 'chains/spectrum_FERMI.dat',
     &   form='formatted',
     &       status='Replace')
!c      msteptb =0.2d0   We use the same as the one in the file stepf.dat
        do k= 1,stepnum
          massvec(k) = 10.d0**( hlog10 - msteptb*(k-1))
c           write(*,*) 'mass ', massvec(k)
         hamwimp = massvec(k) 
          
      do i =1, 42
        ppin= FEmid(i)
           if( i.eq.1)  write(*,*) 'dsepflux ',dsepflux(FEmin(i))
c        FspcV(k,i)=  dsepflux(FEmin(i))
cccccccccccccc UNCOMMENT the following two lines
        call dsfun_intb(dsepflux,FEmin(i),FEmax(i),feps,prec,fluxtemp)
        FspcV(k,i)= fluxtemp / ( FEmax(i) - FEmin(i) )

ccc      to compare to the no integration case
c       fluxtemp = dsepdphidphaaxi(ppin,phiin,how)
c             if( k.eq.2)   write(*,*) FspcV(k,i), fluxtemp
c             if( k.eq.6)   write(*,*) FspcV(k,i), fluxtemp
   
c           write(*,*) 'FspcV ',FspcV(k,i) 
       enddo
      enddo
       DO k = 1 , stepnum
       WRITE(4,44) (FspcV(k,i),i=1,42)
       enddo
!c      DO i = 1 , 42
!c       WRITE(*,*)  FEmid(i)
!c       enddo
44    FORMAT(42(1x,e13.8))
      close(4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   tabulate the step size function spctrum at AMS energy 


       open (unit = 122, file = 'chains/spectrum_AMS.dat',
     &   form='formatted',
     &  status='Replace')
        do k= 1,stepnum
         hamwimp = massvec(k) 
         do i =1, amsnum
        ppin= AMSE(i)
c        AMSspcV(k,i)= dsepflux(AMSE(i))
cccccccccccccc UNCOMMENT the following three lines
        call dsfun_intb(dsepflux,AMSEmin(i),AMSEmax(i),
     &    feps,prec,fluxtemp)
        AMSspcV(k,i)= fluxtemp / ( AMSEmax(i) - AMSEmin(i) )
        enddo
       enddo
       DO k = 1 , stepnum
       WRITE(122,55) (AMSspcV(k,i),i=1,65)
       enddo
55    FORMAT(65(1x,e13.8))
      close(122)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   tabulate the step size function spctrum at pamela energy 


       open (unit = 123, file = 'chains/spectrum_pamela.dat',
     &   form='formatted',
     &  status='Replace')
        do k= 1,stepnum
         hamwimp = massvec(k) 
         do i =1, pamelanum
        ppin= pamelaE(i)
c        pamelaspcV(k,i)= dsepflux( ppin)
cccccccccccccc UNCOMMENT the following three lines
        call dsfun_intb(dsepflux,Rmin(i),Rmax(i),
     &    feps,prec,fluxtemp)
        pamelaspcV(k,i)= fluxtemp / ( Rmax(i) - Rmin(i) )
        enddo
       enddo
       DO k = 1 , stepnum
       WRITE(123,1005) (AMSspcV(k,i),i=1,39)
       enddo
1005    FORMAT(39(1x,e13.8))
      close(123)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc   tabulate the step size function spctrum at hess energy 


       open (unit = 128, file = 'chains/spectrum_hess.dat',
     &   form='formatted',
     &  status='Replace')
        do k= 1,stepnum
         hamwimp = massvec(k) 
         do i =1, hessnum
        ppin= hessE(i)
        hessspcV(k,i)= dsepflux( ppin)
cccccccccccccc UNCOMMENT the following three lines
c        call dsfun_intb(dsepflux,Rmin(i),Rmax(i),
c     &    feps,prec,fluxtemp)
c        hessspcV(k,i)= fluxtemp / ( Rmax(i) - Rmin(i) )
        enddo
       enddo
       DO k = 1 , stepnum
       WRITE(128,1006) (hessspcV(k,i),i=1,hessnum)
       enddo
1006    FORMAT(9(1x,e13.8))
      close(128)













ccccccccccccc  using Fermi Set MAX amplitude


          do k = 1, sdim1
          amptemp1=10.d10
          amptemp2=10.d20
           i=6
                 do j = 6 , 42
         fluxtemp= Fflux(j)+ 4.*Ferrp(j) -bg_electron(j)-bg_positron(j)
          if ( Fspcv(k,j) > 0.)
     &    amptemp2 =   fluxtemp /(2.* FspcV(k,j)  )
             if ( amptemp2< amptemp1) then 
            amptemp1 = amptemp2
             i = j
             end if
                enddo
c            if (amptemp1 < 50.) amptemp1=50.
            stepampmax(k) = amptemp1
            write(*,*) 'step function, max amplitude',i, stepampmax(k)
         end do


ccccccccccccc  using AMS Set MAX amplitude


          do k = 1, sdim1
          amptemp1=10.d10
          amptemp2=10.d20
           i=34
                 do j = 34 ,65 
                 ratiob =( AMSpr(j) + 3.* AMSerr(j))  
         fluxtemp= (ratiob* ( bg_e_AMS(j) + bg_p_AMS(j) ) - bg_p_AMS(j))
     &             / ( 1. - 2.* ratiob);
          if ( AMSspcv(k,j) > 0.)
     &    amptemp2 =   fluxtemp /( AMSspcV(k,j)  )
             if ( amptemp2< amptemp1) then 
            amptemp1 = amptemp2
             i = j
             end if
                enddo
c            if (amptemp1 < 50.) amptemp1=50.
            if (amptemp1 < stepampmax(k))  stepampmax(k) = amptemp1
            write(*,*) 'step function, max amplitude',i, stepampmax(k)
         end do




 1000 format(10(1x,e14.8))
 1001 format(10(1x,e10.4))





         !setting priors
      do  i = 1, sdim1
         spriorran(i,1)=0.
         spriorran(i,2)= stepampmax(i)
      enddo
            
         spriorran(sdim1+1,1)= -0.05d0
         spriorran(sdim1+1,2)= 0.05d0 
 
         spriorran(sdim1+2,1)= -0.2d0
         spriorran(sdim1+2,2)= 0.2d0 

         !no parameters to wrap around
         nest_pWrap=0

 
 
         call spline(FEmid,Fflux,flnum,0.d0,0.d0,Ffluxy2)
         call spline(FEmid,Ferrp,flnum,0.d0,0.d0,Ferry2)
         call nest_Sample


c    read epchi2-sstats.dat file, and calcualte chi2 for the best fit point

      open (unit = 116, file = 
     &'chains/epchi2-stats.dat',
     &   form='formatted')
      open (unit = 117, file = 
     &'chains/dnde_mod1.dat',
     &   form='formatted',status='replace')
      open (unit = 118, file = 
     &'chains/dnde_best.dat',
     &   form='formatted',status='replace')
        do while (tempchar .ne. 'Dim') 
         read( 116,"(A3)") tempchar
         enddo

         write( *,*) 'amplitude and sigma'
        do i = 1, sdim
         read( 116, *) number1,ampbest(i) , res_sigma(i)
         write( 117, *) number1,ampbest(i), res_sigma(i)
         write( *, *) number1,ampbest(i) ,res_sigma(i)
        enddo


        do while (tempchar .ne. 'Max') 
         read( 116,"(A3)") tempchar
c         write( *,*) tempchar
         enddo
         read( 116,"(A3)") tempchar
         write( *,*) 'best fit model'
        do i = 1, sdim
         read( 116, *) number1,ampbest(i) !, res_sigma(i)
         write( 118, *) number1,ampbest(i) !,res_sigma(i)
         write( *, *) number1,ampbest(i) !,res_sigma(i)
         enddo
        close(116)
        close(117)
        close(118)
          amptemp(1:sdim) = ampbest(1:sdim)
           write(*,*) 'chi2 = ', chi2Function(ampbest)
            
       
        end 





       
         real*8   function fermishiftchi2( )
        use params
        implicit none
         real*8 delta,temp1,chi2,slp
         real*8 FEshift(42)
         real*8 Ffluxshift(42)
         real*8 Ferrshift(42)
       integer flnum
       parameter ( flnum = 42)
         integer i,j
              delta  = amptemp(sdim1+1)
            slp = amptemp(sdim1+2)
           do  i = 1 , 42
            FEshift(i) = FEmid(i)* ( 1.+delta)
c         call spline(FEmid,Ferrp,flnum,0.d0,0.d0,Ferry2)
         call splint(FEmid,Fflux,Ffluxy2,flnum,
     &  FEshift(i),Ffluxshift(i))
         call splint(FEmid,Ferrp,Ferry2,flnum,
     &  FEshift(i),Ferrshift(i))
           enddo

          chi2 =0.
           do  i =  42,1,-1
          temp1 =0
              if ( FEmid(i) <20.d0) EXIT
            do  j =1, sdim1
              temp1 = temp1 + amptemp(j)* 2.*FspcV(j,i)
             enddo 
             chi2 = chi2 + (( temp1 +bg_electron(i)
     &  *(( FEmid(i)/10.d0) **slp)  +bg_positron(i)- 
     &   Ffluxshift(i))/Ferrshift(i) )**2

           enddo
          fermishiftchi2 = chi2
            end function


        real*8 function amschi2() 
         use params
          real*8 temp1,chi2,ratio,slp
          integer i,j
            slp = amptemp(sdim1+2)

             chi2=0. 
        do i = 65,1,-1
         temp1 =0.  ! temp1 represents the total flux
          if ( AMSE(i) <20.d0) EXIT
            do  j =1, sdim1
             temp1 = temp1 + amptemp(j)* AMSspcV(j,i)
             enddo
            ratio = ( temp1 + bg_p_AMS(i)) / ( 2.*temp1 +    
     &        bg_e_AMS(i)*((AMSE(i)/10.d0)**slp )+ bg_p_AMS(i)  )
        chi2 = chi2 +( (ratio- AMSpr(i))/ AMSerr(i) )**2
c             write(*,*) i, ratio, AMSpr(i), AMSerr(i),chi2
        enddo
           amschi2 = chi2

        end function
      


        real*8 function pamelachi2() 
         use params
          real*8 temp1,chi2,ratio,slp
          integer i,j
            slp = amptemp(sdim1+2)

             chi2=0. 
        do i = 39,1,-1
         temp1 =0.  ! temp1 represents the total flux
          if ( pamelaE(i) <20.d0) EXIT
            do  j =1, sdim1
             temp1 = temp1 + amptemp(j)* pamelaspcV(j,i)
             enddo
         chi2 = chi2 + (( temp1 +    
     &      bg_e_pamela(i)* ( pamelaE(i)/10.d0)**slp- 
     &     pamelaflux(i))/pamelaerr(i) )**2

        enddo
           pamelachi2 = chi2

        end function



ccccccccccccccc HESS  chi2 function 
        real*8 function hesschi2() 
         use params
          real*8 temp1,chi2,slp,fluxtemp1
          integer i,j
            slp = amptemp(sdim1+2)

             chi2=0. 
        do i = 9,1,-1
         temp1 =0.  ! temp1 represents the total flux
          if ( hessE(i) <20.d0) EXIT
            do  j =1, sdim1
             temp1 = temp1 + amptemp(j)* hessspcV(j,i)
             enddo
             fluxtemp1 = temp1 + bg_e_hess(i)* 
     &   ( hessE(i)/10.d0)**slp + bg_p_hess(i)
            if ( fluxtemp1 > hessflux(i)) then
         chi2 = chi2 + ((fluxtemp1 - hessflux(i))/hesserr(i))**2
             end if

        enddo
           pamelachi2 = chi2

        end function






 

       real*8 function fermishift(delta)
         implicit none
            real*8  fermishiftchi2
           real*8  delta
            real*8 chi2
          real*8 sigma
           sigma = 0.08

              chi2 = fermishiftchi2(delta)
           fermishift=chi2* dexp( - delta*delta/ ( 2.*sigma*sigma))/
     &    (dsqrt(2.d0* 3.1415926)* sigma)
           
        end function
              
          
         real*8   function int_fermishift(a,b)
           implicit none
            real*8  fermishift, a,b,eps,prec
            real*8 res
             
             eps =1.d-3
             prec =1.d-3
c              write( *,*) 'fermishift a ',fermishift(a)
c              write( *,*) 'fermishift b ',fermishift(b)
c              write( *,*) 'fermishift 0 ',fermishift(0.d0)
c              write( *,*) 'fermishift 0.3 ',fermishift(0.3d0)
                res  = fermishift(0.d0)
          call dsfun_int( fermishift,a,b,eps, prec, res)
             int_fermishift =res
         end function
            

      




               include 'funs_ep.f'
cc=====================================================================

