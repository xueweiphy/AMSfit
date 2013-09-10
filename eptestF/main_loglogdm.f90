      include 'user-dummy-eptest.h'
      include 'user-not-dummy-eptest.h'

      program main
      use nestwrapper
      use params
      use HydroDensity
      use InterpolData
      use DataChi2
         

ccc routine to find the DM spectrum to fit Fermi, HESS, AMS, PAMELA
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
      real*8  engvecmin(121)
      common / engveccom/ massvec,spcvec,engvec, engvecmin
    
      real*8 dsepflux,dsf_int
ccc
ccc for the user-defined ep initial spectrum
ccc
c      integer modeep
c      common/dsepdndpiusercom/modeep

      real*8 mstep1,engup1,alphaloglog
      common /massstep/ mstep1,engup1,alphaloglog
      real*8 listalphacom(stepnum)
      common listalphacom

ccc
      integer nbinsep
      real*8 deltalogS(100),deltalogp(100),logSendpoint,logpendpoint
      common/dsepdndpiuser2com/deltalogS,deltalogp,logSendpoint,
     &  logpendpoint,nbinsep


ccc      step function amplitude
      real*8 stepamp(sdim1),stepampmax(sdim1)
      real*8  amptemp1, amptemp2
      real*8 ampbest(sdim), res_sigma(sdim),cube(sdim)
      integer number1
      character(len=3) tempchar
c          common/stepcom/ stepamp,stepampmax,amptemp1,amptemp2,
c     &    ampbest,res_sigma,cube,number1,tempchar 

ccc
      real*8 alphap,emaxp,eminp,ampdnde,normeng
      common /dsepdndpiuser4com/alphap,emaxp,eminp,ampdnde,normeng
      real*8 dsepdndpi_user
      external dsepdndpi_user
      real*8 pmin,pmax,eps,res,par,dsdndepulsar,age,dist,etot,norm


      integer nbgmodel
      real*8 Rsun
      common/dsnpaxiusercom/Rsun,nbgmodel



      real*8 xmin,xmax,deltav,result,dsepgreenaxi,dsepgreenaxitab,pp,
     &  dsepdndpaxi

      real*8 normprim   ! primary electron normaliztion
      common/normprimcom/normprim
      real*8 normprim0, alphap0
      integer flagnorm
      
      integer intflag

      
      character*12 spcnamelabel
      integer binnumle
      real*8   logenglist(50),englist(50),amplitudeR(50)
      common    /englistcom/ binnumle,logenglist,amplitudeR
       
      real*8 normlist_test(sdim1),alphalog0(sdimlog)
      real*8 flux_F_test(flnum),flux_F_test2(flnum)
      real*8 flux_A_test(amsnum), flux_A_test0(amsnum)
      real*8 flux_AMSep_test(AMSeplnum)
      integer flagoutr
      real*8 ngas
      real*8 engtest 


      call dsinit


      call AMSDataIni
      call readdatafiles
      
c      write(*,*) 'AMS ep data'
c      do i = 1, AMSeplnum
c         write(*,*) AMSepE(i),AMSepF(i), AMSepErru(i), AMSepErrd(i)
c      enddo


      call nTotal_gal(8.5d0, 0.d0, ngas)
      write(*,*) 'ngas',ngas
     
      call SecSourceIni()
 
      do i  = 1, 15
         write(*,*) 
     &  AMSelectronE(i-15+ AMSelnum),
     &  AMSelectronF(i-15+ AMSelnum),
     &  AMSelectronErru(i-15+ AMSelnum),
     &  AMSelectronErrd(i-15+ AMSelnum)
      enddo
      do i = 1, 15
         engtest = i*1.
         write(*,*) engtest, AMSProton_Flux( engtest),
     &     SecPositronSource(engtest, 8.5d0,0.d0)
      enddo

      binnumle  = stepnum
      intflag =0

      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
ccc    recalculate the primary and seconard electron and poisitron
ccc

ccc halo model -> replaced by background source distribution
ccc
      npaxidefault=.false.  ! do not link to spherical distribution
ccc
ccc energy loss model M1 in Table 2 of Delahaye et al. A&A 524, A51 (2010)
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


      diffrcep=0.001d0    ! inner radial cut in ep routines in kpc, this value
                          ! should not matter
      epcraxitag='dskra'
      call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep) ! add label for 
                          ! book keeping and correct link to tabulated green
                          ! function

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

ccc
ccc primary electrons, Ferriere radial distribution
ccc
      nbgmodel=2
      Rsun=8.5d0
      nptag='ferriere'   ! tag needed for tabulation of the green function
      nnptag=8



ccc check the green function:
ccc

      how=4
      deltav=1.d0
      result=dsepgreenaxitab(Rsun,0.d0,DeltaV,how)
      write(*,*) deltav,result



ccc
ccc to get equilibria spectra, namely the integral of the green function 
ccc on the choose spectrum you have to fool the code and let it think
ccc it is some WIMP of given mass and specify the "line" branching ratio
ccc
      hamwimp=1.d4 ! this is the maximum compatible with vvarnum tabulation
                   ! if you need it larger than this change the tabulation
      habr(15)=0.d0



ccc NOTE: NOW THERE IS A FAKE ENERGY SPECTRUM!!!!!!!!!!!!!!!!!!!!!
ccc the one introduced for pulsars for modeep=5; FIX THIS
ccc

ccc
ccc dn/dE is power law + cutoff in energy
ccc dn/dp= dn/dE * dE/dp
ccc the normalization is such that int_Emin^infinity dE E dn/dE = 1
ccc hence result should be scaled up by the total energy into positrons
ccc
      epspecdef=.false.
      modeep=9
      alphap=2.7d0
      emaxp=30000.d0
      ampdnde = 0.004/( 3.d8/4./3.141593)
      normeng = 33.d0
c     eminp=mel
      eminp=0.1d0 ! ????? this is what they quote for SNR not pulsar


ccc      set initial value of normprim
      normprim = 5.03d-11


      xmin=1.d0
      xmax=1000.d0

c      do i=0,10
c        pp=dexp(dlog(xmin)+(dlog(xmax)-dlog(xmin))/10.d0*i)
c        result=dsepdndpaxi(Rsun,0.d0,pp,ivopt,how)
c        write(*,*) i,pp,result,result*pp**2,result*pp**3
c      enddo




ccc
ccc to get the flux fix the conversion from number density to flux,
ccc adding proper units, normalizations, so on so forth !!!!!
ccc NOTE: here it is not possible to use dsepdphidphaaxi, since it really
ccc refers to WIMPS
ccc


ccc        Use PAMELA data at 30.6 Gev to normalize the flux


      alphap0= 2.7d0       
      normprim0 =1.d0
      flagnorm =1  ! 1 is to reset the normalization
      do i = 1, stepalpha
         alphap0 = 1.9 + 0.1*real(i)
         prim_alpha(i) = alphap0
         call setprie_exp( alphap0,normprim0,flagnorm)
         prim_F_table(i,:) =  bg_prie_F
         prim_A_table(i,:) =  bg_prie_A
         prim_P_table(i,:) =  bg_prie_P
         prim_H_table(i,:) =  bg_prie_H
         
         Prim_AMSelectronAlpha(i,:)= Prim_AMSelectron
         Prim_AMSpositronAlpha(i,:)= Prim_AMSpositron
         Prim_AMSepAlpha(i,:)= Prim_AMSep
      enddo
      call SplinePrimAMS

      call spline_prim
      alphap0 = 2.70
      call splint_prim(alphap0,normprim0)
      call SplintPrimAMS(alphap0,normprim0)
      do i = 1, 40
         write(*,*)i,FEmid(i),bg_prie_F(i), AMSelectronE(i),
     &      Prim_AMSelectron(i)
      enddo

      call setsecep_exp   ! calculate the secondary electron and positr    on 
      call addprisece



   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc DM halo model:
ccc
      npaxidefault=.true.  ! link to spherical distribution
      npsphdefault=.true.  ! link to rho^2
      rhodefault=.false.   ! do not link to the default distribution but to
                           ! but to dsnpsph_user

ccc
ccc choose a halo profile and make a dummy call to it set the 
ccc galactocentric distance and the inner truncation radius which are
ccc hardcoded in dsrhosph_user but can be used as global variables
ccc
      
      nptag='pbnfw'   ! this is the reference NFW profile used in the
                      ! pbar paper, see user-not-dummy-ep/dsrhosph_user.f
      nnptag=5
       
      dummyres=dsrhosph_user(1.d0)
      epgraxiloadck=.true.
       write(*,*) 'reload tabulation'
      epcraxitag='dskra_test'
      call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep) 
       
      epgraxiloadck=.false.


ccc
ccc choose the diffusion coefficient and halo height 
ccc

c      include 'propmod/pbkra.f'  ! reference kraichan model in the pbar paper
      diffrcep=0.001d0    ! inner radial cut in ep routines in kpc, this value
                          ! should not matter


ccc
ccc model "H" for energy loss parameters in the mean field model
ccc
c      ivopt=1  ! the matching between  pp and the variable v is done 
               ! assuming fully general momentum loss rate and spatial
               ! diffusion coefficient, with a numerical integral + a
               ! tabulation involved
ccc for ivopt=1:
ccc
ccc some mean from moskalenko talk or porter paper
ccc
ccc      Here I use the same star light, comment the follwinglin out
c      starlightmean=1.d0-0.25d0 ! optical + IR, eV cm^-3
ccc
ccc some 'local' mean in 8*dexp(-(r-r0)/50kpc) *dexp(-|z|/3kpc)
ccc random + regular, strong talk 
ccc
c      bfieldmean=6.d0      ! \muG

ccc
ccc diffusion time for a few values of the pbar kinetic energy:
ccc
      how=4       ! green function read from table if it exists
      phiin=0.5d0   ! forced field solar modulation parameter 

ccc
ccc
ccc step function for the input spectrum:
ccc
      epspecdef=.false.
ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccc       change modeep 
ccccccccc   3--- step function
cccccccc    4--- theta function  
ccc   DO NOT Set modeep to 8 here, use 3 or 4 to set the scan range
      spcnamelabel = 'dmstep'
      modeep= 3
      hasv=3.d-26 ! cm^3 s^-1 -> you still need this because dsepdphidphaaxi
                  ! is taken to be proportional to hasv/hamwimp**2
c      phiin=0.d0   ! look at the interstellar

c       habr(15)=0.d0  ! this sets the monochromatic term at E=hamwimp,
                     ! for habr(15)=1.d0 the spectrum is just bypassed
       


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc tabulate the step size function' spectrum.
cccc  mass starting from few TeV to 
cccc  the mass smaller 10 GeV
ccc   energy from  1 GeV to 1 TeV
ccc

      do k =1, stepnum
         engvec(k) = 10.d0**( hlog10- engsteptb*(k-1))
         engvecmin(k) = 10.d0**( hlog10- engsteptb*(k))
         amplitudeR(k) = 0.d0
         logenglist(stepnum+1-k)= hlog10- engsteptb*(k-1)
      enddo



      hamwimp = 1.d4

      mstep1 = msteptb
      do k =1,  sdim1 
         massvec(k) = 10.d0**( hlog10- msteptb*(k-1))
         hamwimp = massvec(k) 
         engup1 = massvec(k) 
         do i =1, 146
            ppin= kra_eng(i)
            if ( ppin > massvec(k) ) then
               spcvec(k,i) =0.
            else 
               spcvec(k,i) = dsepdphidphaaxi(ppin,phiin,how)
            end if
         enddo 
      enddo

c      call writekra

      feps =1.d-3
      prec = 1.d-2
      call spectrumdm_data(intflag,0)
      write(*,*) 'modedp = ', modeep





      
            

ccccccccccccc  using Fermi Set MAX amplitude
        do k = 1, sdim1
           amptemp1=10.d10
           amptemp2=10.d20
           i=6
           do j = 6 , 42
              fluxtemp= Fflux(j)+ 10.*Ferrp(j) -
     &                  bg_electron(j)-bg_positron(j)
              if ( Fspcv(k,j) > 0.)
     &           amptemp2 =   fluxtemp /(2.* FspcV(k,j)  )
              if ( amptemp2< amptemp1) then 
                 amptemp1 = amptemp2
                 i = j
              end if
           enddo
c          if (amptemp1 < 50.) amptemp1=50.
           stepampmax(k) = amptemp1
           write(*,*) 'step function, max amplitude',i, stepampmax(k)
        end do


ccccccccccccc  using AMS Set MAX amplitude


       do k = 1, sdim1
          amptemp1=10.d10
          amptemp2=10.d20
          i=34
          do j = 34 ,65 
             ratiob =( AMSpr(j) + 10.* AMSerr(j))  
             fluxtemp= (ratiob* ( bg_e_AMS(j) + bg_p_AMS(j) )
     &             - bg_p_AMS(j))  / ( 1. - 2.* ratiob);
c             if (k .eq. 2)
c     &          write(*,*) j, bg_e_AMS(j), bg_p_AMS(j), fluxtemp 
             

             if ( AMSspcv(k,j) > 0.)
     &          amptemp2 =   fluxtemp /( AMSspcV(k,j)  )
             if ( amptemp2< amptemp1) then 
                amptemp1 = amptemp2
                i = j
             end if
          enddo
c            if (amptemp1 < 50.) amptemp1=50.
          if (amptemp1 < stepampmax(k))  stepampmax(k) = amptemp1
          write(*,*) 'step function, max amplitude',i,k, stepampmax(k)
       end do


      write(*,*) " FERMI energy list"
      write(*,*) FEmid 


 1000 format(10(1x,e14.8))
 1001 format(10(1x,e10.4))

c      epgraxiloadck=.true.
c      npaxidefault=.false.
c      epspecdef=.false.
c      hamwimp = 1.d4
c      epcraxitag='dskra'
c      call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep) 
c      write(*,*) 'for test1'
c      call setprie_exp( alphap0,normprim0,flagnorm)
c      write(*,*) 'alpha = ',alphap0
c      epgraxiloadck=.false.
c      call setsecep_exp
c      call addprisece 
c      write(*,*) 'modeep = ',modeep

      hamwimp = 1.d4

      normlist_test= stepampmax  
      normlist_test(1)=   23.506706724698837
      normlist_test(2)=    2.8557846343755360
      normlist_test(3)=   59.935765492490127
      normlist_test(4)=    168.72288432662558
      normlist_test(5)=   564.21479693851961
      normlist_test(6)=   2095.4803217067756
      normlist_test(7)=  7051.4012268872448
c      normlist_test(8)=   23283.657433871729
c      normlist_test(9)=   114762.27063383316

      modeep = 10
      amptemp(1:sdim1) = normlist_test
      call spectrumdm_data(intflag,0)
      do  i = 1, flnum
         flux_F_test2(i) = FspcV(1,i)
      enddo

      do  i = 1, amsnum
         flux_A_test0(i) =AMSspcV(1,i)
      enddo
         
      write(*,*) "hamwimp",hamwimp
      
c      flux_F_test2 = FspcV(1,:)
c      flux_A_test0 = AMSspcV(1,:)

      modeep = 11
      mstep1 = msteptb
         do i = 1, 61
            alphaloglog  = (real(i)-31.0)*0.1 /mstep1
            alphaloglist(i) = alphaloglog
            listalphacom = alphaloglog
            call spectrumdm_data(intflag,0)
            dm_A_table(:,:,i)=AMSspcV(1:sdimlog,:)
            dm_F_table(:,:,i)=FspcV(1:sdimlog,:)
            dm_P_table(:,:,i)=pamelaspcV(1:sdimlog,:)
            dm_H_table(:,:,i)=hessspcV(1:sdimlog,:)
            dmps_AMSpositronAlpha(:,:,i)=AMSpositronSpcV(1:sdimlog,:)
            dmps_AMSelectronAlpha(:,:,i)=AMSelectronSpcV(1:sdimlog,:)
            dmps_AMSepAlpha(:,:,i)=AMSepSpcV(1:sdimlog,:)
         enddo
     
 
c      write(*,*) 'test AMSepAlpha(:,:,:)'
c      do i = 1, 40
c         write(*,*) i, FEmid(i), dm_F_table(2,i,4), AMSepE(i),
c     &      dmps_AMSepAlpha(2,i,4)
c      enddo

c      write(*,*) dmps_AMSepAlpha(2,:,4) 
c      write(*,*) dmps_AMSepAlpha(2,:,4) 

      call spline_spcloglog! For Fermi, AMSper, PAMELA and HESS data
      call SplineDmpsAMS   ! For AMS data
       
     
      call alpha_loglog(normlist_test,alphalog0,flagoutr)
c      write(*,*) 'alphalog' 
c      write(*,*) alphalog0 
c      write(*,*) 'normlist' 
c      write(*,*) normlist_test 
     
      call splint_spcloglog(normlist_test,alphalog0)
      call splintDmpsAMS(alphalog0)
      flux_F_test =0.d0
      flux_AMSep_test =0.d0
      do i = 1 , sdimlog
        flux_F_test = flux_F_test + normlist_test(i)*FspcV(i,:)
      enddo
      flux_A_test =0.d0
      do i = 1 , sdimlog
        flux_A_test = flux_A_test + normlist_test(i)*AMSspcV(i,:)
         flux_AMSep_test = flux_AMSep_test+ normlist_test(i)*
     &      AMSepSpcV (i,:)
      enddo
       
c      do i = 1, 40
c         write(*,*) i, FEmid(i), flux_F_test(i), AMSepE(i),
c     &      flux_AMSep_test(i)
c      enddo

      write (*,*) chi2AMSep() 

c      return
      

      open (unit = 102, file = 'Root/spectrum_F_test.dat',
     &  form='formatted' , status='Replace')

      open (unit = 103, file = 'Root/spectrum_A_test.dat',
     &  form='formatted' , status='Replace')

      write(*,*) ' write spectrum to spectrum_stepf'

         
      write(103,1030) AMSE 
      write(103,1030) flux_A_test0
      write(103,1030) flux_A_test

      write(102,1020) FEmid 
      write(102,1020) flux_F_test2
      write(102,1020) flux_F_test
1020  FORMAT(42(1x,e14.8) )
1030  FORMAT(65(1x,e14.8) )

      listalphacom(:sdimlog) = alphalog0
      listalphacom(sdimlog+1) = 0.d0
      call spectrumdm_data(intflag,0)
      flux_F_test =0.d0
      do i = 1 , sdimlog
        flux_F_test = flux_F_test + normlist_test(i)*FspcV(i,:)
      enddo
      write(102,1020) flux_F_test


       close(102)
      close(103)
c      write(*,*) FspcV(5,:) 
c      write(*,*) AMSspcV(1,:) 
c      write(*,*) AMSspcV(2,:) 
c      write(*,*) hessspcV(5,:) 
c      modeep =11
c      mstep1 =  msteptb*(sdim1-1)  
c      alphaloglog =0.1
c         hamwimp =  10.d0 ** hlog10 
c         write(*,*) 'hamwimp',hamwimp
c         engup1 =  hamwimp 
c      do i =70, 90
c         ppin= kra_eng(i)
c         spcvec(1,i) = dsepdphidphaaxi(ppin,phiin,how)
c         write(*,*) i, dsepdndpi_user(ppin)*((10000./hamwimp)**2)
c         write(*,*) i, ppin,dsepdndpi_user(ppin)
c      enddo
c         ppin  = 10.** logenglist(sdim1-1) 
c         write(*,*) ppin, dsepdndpi_user(ppin)
c      write(*,*) spcvec(1,70:90)*((10000./hamwimp)**2)
c      write(*,*) spcvec(1,70:90)
c      modeep = 8
c      do i = 1, sdim1
c         amplitudeR(sdim1-i+1)=10.**( (real(i) -1.)*0.1*msteptb)
c      enddo
c      do i =70, 90
c         ppin= kra_eng(i)
c         spcvec(1,i) = dsepdphidphaaxi(ppin,phiin,how)
c         write(*,*) i, ppin,dsepdndpi_user(ppin)
c      enddo
c      ppin  = 10.** logenglist(sdim1-1) 
c      write(*,*)  ppin,dsepdndpi_user(ppin),amplitudeR(sdim1-1)
c      write(*,*) spcvec(1,70:90)





         !setting priors
   
        do i = 1, sdim1
           spriorran(i,1)= 0.
           spriorran(i,2)= stepampmax(i) *100.d0
        enddo
            
        spriorran(sdim1+1,1)= -0.05d0 !shift
        spriorran(sdim1+1,2)= 0.05d0 
 
        spriorran(sdim1+2,1)= 2.2d0   !primary injection
        spriorran(sdim1+2,2)= 3.3d0 

        spriorran(sdim1+3,1)= 0.1d0   !primary norm
        spriorran(sdim1+3,2)= 10.d0 

        !no parameters to wrap around
        nest_pWrap=0

 
 
         call spline(FEmid,Fflux,flnum,0.d0,0.d0,Ffluxy2)
         call spline(FEmid,Ferrp,flnum,0.d0,0.d0,Ferry2)
         
         
         modeep =11
         hamwimp = 1.d4

         call nest_Sample


c    read epchi2-sstats.dat file, and calcualte chi2 for the best fit point

         open (unit = 116, file = 
     &      'chains/epchi2-stats.dat',
     &      form='formatted')
         open (unit = 117, file = 
     &      trim(prefix_f)//'dnde_mod1.dat',
     &      form='formatted',status='replace')
         open (unit = 118, file = 
     &      trim(prefix_f)//'dnde_best.dat',
     &      form='formatted',status='replace')
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
c         amptemp(1:sdim) = ampbest(1:sdim)
         amptemp = ampbest
         alphap0 =  ampbest(sdim -1)
         normprim0 = ampbest(sdim)
         flagnorm =1
c         call setprie_exp( alphap0,normprim0,flagnorm)
         call splint_prim( alphap0,normprim0)
         call SplintPrimAMS(alphap0,normprim0)
         call setsecep_exp   ! calculate the secondary electron and positron 
         call addprisece 
         intflag =1
         if ( modeep .eq. 10) then
            call spectrumdm_data(intflag,1)
         else if ( modeep .eq. 11) then
            call alpha_loglog(ampbest(1:sdim1),alphalog0,flagoutr)
            call splint_spcloglog(ampbest(1:sdim1),alphalog0)
            call SplintDmpsAMS(alphalog0)
         endif

         write(*,*) 'chi2 = ', chi2Function(ampbest)


      hamwimp = 1.d4
      mstep1 = msteptb
      spcvec =0 
      do k =1,  sdimlog 
         massvec(k) = 10.d0**( hlog10- msteptb*(k-1))
         hamwimp = massvec(k) 
         engup1 = massvec(k) 
         alphaloglog = alphalog0(k)
         listalphacom = alphaloglog
         do i =1, 146
            ppin= kra_eng(i)
            if ( ppin  .lt. hamwimp)  then
                spcvec(k,i) = dsepdphidphaaxi(ppin,phiin,how)
            else
                spcvec(k,i) = 0.d0
             endif
         enddo
      enddo

      write(*,*) "prior "
      write(*,*)  spriorran(:,2)

      call writekra
      call writedata
      call   priep_kra(alphap0,normprim0, flagnorm)




            
       
       end 




