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


ccccc   from pulsar test file below
         real*8 res
        real*8 dsepdphidp_setpul,deltat
        character*100 popfile
      common/pulsarfluxc/ popfile,deltat
       

      real*8 alphap,emaxp,eminp
      common /dsepdndpiuser4com/alphap,emaxp,eminp,ampdnde,normeng
      external dsepdndpi_user


      real*8 parres(100)
      integer parnum(100)
      common/pulsarsetcom/parres,parnum
ccc
      real*8 dist(100),age(100),edot(100),Weplus(100)
      integer num
      character*10 name1(100),name2(100)
      common/pulsarpopcom/dist,age,edot,Weplus,num,name1,name2
ccc
      real*8 mel
      parameter (mel=0.000510999907d0) !e-mass in GeV
ccc
      character*10 name1d,name2d
      character*100 scr
      real*8 distd,aged,etotd,fraction,dsepdphidp_1plpul,taudec,vexp
      character*7 dummy
      integer ii,jj
ccccc   from pulsar test file above


 



      real*8 radius,res1,res2,dsrhosph_user

      real*8 dummyres

      real*8 mchiin,sigmavin,ppup,pplow,resmat(4,100),dsepdndpi
      integer k,chin,chinvec(4),nbins
      real*8 mstepvec(4),resvec(4,100),engstep
      real*8 massvec(31),spcvec(31,146),engvec(121)
      real*8  engvecmin(121)
       common / engveccom/ massvec,spcvec,engvec, engvecmin
    
      real*8 dsepflux,dsf_int
      real*8 dsepfluxps
ccc
ccc for the user-defined ep initial spectrum
ccc
      integer modeep
      common/dsepdndpiusercom/modeep

      real*8 mstep1,engup1
      common /massstep/ mstep1, engup1

ccc
      integer nbinsep
      real*8 deltalogS(100),deltalogp(100),logSendpoint,logpendpoint
      common/dsepdndpiuser2com/deltalogS,deltalogp,logSendpoint,
     &  logpendpoint,nbinsep

ccc
       integer  binnumle
       real*8   logenglist(50),englist(50),amplitudeR(50)
       common    /englistcom/ binnumle,logenglist,amplitudeR



ccc      step function amplitude
          real*8 stepamp(sdim1),stepampmax(sdim1)
          real*8  amptemp1, amptemp2
          real*8 ampbest(sdim), res_sigma(sdim),cube(sdim)
          integer number1
          character(len=3) tempchar
c          common/stepcom/ stepamp,stepampmax,amptemp1,amptemp2,
c     &    ampbest,res_sigma,cube,number1,tempchar 

          integer intflag  ! whether to do integration 

      call dsinit



      call readdatafiles

c       binnumle = sdim1 
       binnumle = stepnum

ccc
ccc energy loss model M1 in Table 2 of Delahaye et al. A&A 524, A51 (2    010)
ccc
      nisrf=6
      t0isrf(1)=2.726d0
      Uradisrf(1)=1.d0
      t0isrf(2)=33.07d0
      Uradisrf(2)=4.5d-5
      t0isrf(3)=313.32d0
      Uradisrf(3)=1.2d-9
      t0isrf(4)=3249.3d0
      Uradisrf(4)=7.03d-13
      t0isrf(5)=6150.4d0
      Uradisrf(5)=3.39d-14
      t0isrf(6)=23209.0d0
      Uradisrf(6)=8.67d-17
      bfieldmean=1.d0
ccc
ccc match in case you want to take only the Thomson regime
ccc again model M1 in Table 2
ccc
      starlightmean=(25.4d0+5.47d0+37.d0+22.9d0+11.89d0)*1.d-2
ccc
ccc "med" propagation model in the Annecy/Torino group
ccc
      kdiffdef=.true. ! use the default ds hardcoded version of the 
                      ! diffusion coefficient, dskdiff_def
      nkdiff=1
      kdiffrig0=1.d0
      k0halo=0.0112d0 ! kpc^2/Myr
      k0halo=k0halo*(3.0856d21)**2/(1.d6*31556926d0)*1.d-27 ! 10^27 cm^2 s^-1
      k0gasdisc=k0halo
      write(*,*) 'K0 : ',k0halo
      kdiffdeltalow=0.d0  ! spectral index below the break - dummy variable
      kdiffdelta=0.7d0    ! spectral index above the break
      kdiffeta=1.d0    ! multiply the diffusion coefficient by beta**kdiffeta
ccc
ccc other propagation parameters (pb, db & ep):
ccc
      diffhh=4.d0     ! 1/2 of vertical size of the diff. region in kpc
ccc
ccc other propagation parameters (pb & db only):
ccc
      diffRh=30.d0    ! radial size of the diff. region in kpc
      diffhg=0.1d0    ! half thickness of the disk in kpc
      diffng=1.d0     ! density of hydrogen in the disc in cm^-3
      diffnh=0.d0     ! density of hydrogen in the halo in cm^-3
      diffcvel=0.d0   ! galactic wind velocity in km/s
      pbcraxitag='dsmed'  

      ivopt=1
ccc
ccc Thomson regime:
ccc
c      nisrf=0
ccc
ccc full energy losses:
ccc
      nisrf=6

ccc
ccc this is a check on whether v & u tabulations needs to be reloaded
ccc consider adding this check at each u & v call
ccc
      call vvarnumsetup
      call uvarnumsetup
  
       epspecdef=.false.
       
 


ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccc       change modeep 
ccccccccc   6--- step function
cccccccc    7--- theta function
cccccccc    8--- linear interpolation function
      modeep=6 ! Here you should NOT set modeep to 8, 
               ! since I use 6 or 7 to set scan region

      phiin=0.0d0   ! forced field solar modulation parameter 
      deltat=0.d0

ccc
ccc this is the catalogue for our analysis !!!!!!!!!!!!!!!!!
ccc TAKE IT FROM HERE
ccc
      popfile='eptestF/pulsar-cat130508.txt'
 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc tabulate the step size function. mass starting from few TeV to 
cccc  the mass smaller 10 GeV
ccc   energy from  1 GeV to 1 TeV
ccc
      do k =1, stepnum 
         engvec(k) = 10.d0**( hlog10- engsteptb*(k-1))
         engvecmin(k) = 10.d0**( hlog10- engsteptb*(k))
         amplitudeR(k) = 0.d0
         logenglist(stepnum+1-k)= hlog10- engsteptb*(k-1) 
      enddo

      feps =1.d-3
      prec = 1.d-2
      call spectrumkra
c          call writekra


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc tabulate spectrum . mass from ~1 TeV to ~10 GeV
ccc   energy are at Fermi, AMS, PAMELA, HESS data points
ccc
       intflag = 1
       call spectrumdata( intflag) 
c      call writedata
       

       amptemp(1:stepnum) = 1.d1
       WRITE(*,*)  'mode 7 -----------------'  
            modeep=7
         call spectrumdata(0)
       DO k = 1  ,1
c       WRITE(*,*) (spcvec(k,i),i=1,146)  
c       WRITE(*,*) (spcvec(k,i),i=74,75)  
         do i = 1, 42 
c       WRITE(*,*) AMSE(i), AMSspcV(k,i)  
       WRITE(*,*) FEmid(i), FspcV(k,i)  
           enddo
       enddo
22     FORMAT(146(1x,e13.8))

       WRITE(*,*)  'mode 8 -----------------'  
       amptemp=1.e-30 
            do i = 1, stepnum 
       amptemp(i) =  1.d0 
           enddo
       
 
            modeep=8
         call spectrumdata(0)
c         call spectrumkra
c       WRITE(*,*) (spcvec(k,i),i=1,146)  
c       WRITE(*,*) (spcvec(k,i),i=74,75)  
             k=1
         do i = 1, 42 
        
       WRITE(*,*) FEmid(i), FspcV(1,i)  
c       WRITE(*,*) AMSE(i), AMSspcV(k,i)  
         enddo
c       WRITE(*,*) (spcvec(2,i),i=70,140)  


c          return
            modeep=7
        call spectrumdata( intflag) 
ccccccccccccc  using Fermi Set MAX amplitude


          do k = 1, sdim1
          amptemp1=10.d10
          amptemp2=10.d20
           i=7
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
            write(*,*) 'step function, max amplitude',i, k,stepampmax(k)
         end do

 
ccccccccccccc  using AMS Set MAX amplitude


          do k = 1, sdim1
c          amptemp1=10.d10
c          amptemp2=10.d20
           i=34
                 do j = 34 ,65 
                 ratiob =( AMSpr(j) + 3.* AMSerr(j))  
         fluxtemp= (ratiob* ( bg_e_AMS(j) + bg_p_AMS(j) ) - bg_p_AMS(j))
     &             / ( 1. - 2.* ratiob);
c                write(*,*) j, bg_e_AMS(j) , bg_p_AMS(j)
          if ( AMSspcv(k,j) > 0.)
     &    amptemp2 =   fluxtemp /( AMSspcV(k,j)  )
             if ( amptemp2< amptemp1) then 
            amptemp1 = amptemp2
             i = j
             end if
                enddo
c            if (amptemp1 < 50.) amptemp1=50.
            if (amptemp1 < stepampmax(k))  stepampmax(k) = amptemp1
            write(*,*) 'step function, max amplitude',i,k, stepampmax(k)
         end do




 1000 format(10(1x,e14.8))
 1001 format(10(1x,e10.4))




         !setting priors
      do  i = 1, sdim1
         spriorran(i,1)=0.1d-30
c         spriorran(i,2)= 10000.*stepampmax(i)
         spriorran(i,2)=10.
      enddo
         spriorran(sdim1,2)=700. 
         spriorran(sdim1-1,2)=150. 
         spriorran(sdim1-2,2)=150. 
         spriorran(sdim1-3,2)=150. 
         spriorran(sdim1-4,2)=150. 
            
         spriorran(sdim1+1,1)= -0.05d0
         spriorran(sdim1+1,2)= 0.05d0 
 
         spriorran(sdim1+2,1)= -0.2d0
         spriorran(sdim1+2,2)= 0.2d0 

         !no parameters to wrap around
         nest_pWrap=0

 
 
            
         call spline(FEmid,Fflux,flnum,0.d0,0.d0,Ffluxy2)
         call spline(FEmid,Ferrp,flnum,0.d0,0.d0,Ferry2)




          modeep =8
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

         call spectrumkra
          call writekra
          intflag = 0
        call spectrumdata( intflag) 
        call writedata
              
           write(*,*) 'chi2 = ', chi2Function(ampbest)
            
       
        end 














