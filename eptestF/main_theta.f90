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
      real*8  engvecmin(121)
      common / engveccom/ massvec,spcvec,engvec, engvecmin
    
      real*8 dsepflux,dsf_int
ccc
ccc for the user-defined ep initial spectrum
ccc
      integer modeep
      common/dsepdndpiusercom/modeep

      real*8 mstep1, engup1
      common /massstep/ mstep1,engup1

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

      
      character*12 spcnamelabel



      call dsinit


      call readdatafiles


      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
ccc    recalculate the primary and seconard electron and poisitron
ccc

ccc halo model -> replaced by background source distribution
ccc
      npaxidefault=.false.  ! do not link to spherical distribution
ccc
ccc energy loss model M1 in Table 2 of Delahaye et al. A&A 524, A51 (201    0)
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
      write(*,*) prefix_f//"hhahah"
      return
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
      call setprie_exp( alphap0,normprim0,flagnorm)
      call setsecep_exp   ! calculate the secondary electron and positron 
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
      spcnamelabel = 'dmstep'
      modeep=4
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




      mstep1 = msteptb
c      do k =1, stepnum 
      do k =1,  sdim1 
         massvec(k) = 10.d0**( hlog10- msteptb*(k-1))
c         hamwimp = massvec(k) 
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
      
      call writekra

      feps =1.d-3
      prec = 1.d-2
      call spectrumdm_data(0)
      call writedata







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
          write(*,*) 'step function, max amplitude',i, stepampmax(k)
       end do




 1000 format(10(1x,e14.8))
 1001 format(10(1x,e10.4))

      epgraxiloadck=.true.
      npaxidefault=.false.
      epspecdef=.false.
      hamwimp = 1.d4
      epcraxitag='dskra'
      call dsaddlabelepaxi(epcraxitag,diffhh,diffrcep) 
      write(*,*) 'for test1'
      call setprie_exp( alphap0,normprim0,flagnorm)
      write(*,*) 'alpha = ',alphap0
      epgraxiloadck=.false.
      call setsecep_exp
      call addprisece 
      write(*,*) 'modeep = ',modeep




         !setting priors
   
        do i = 1, sdim1
           spriorran(i,1)= 0.
           spriorran(i,2)= stepampmax(i)
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
         call nest_Sample


c    read epchi2-sstats.dat file, and calcualte chi2 for the best fit point

         open (unit = 116, file = 
     &      'chains/epchi2-stats.dat',
     &      form='formatted')
         open (unit = 117, file = 
     &      'chains/dnde_mod1.dat',
     &      form='formatted',status='replace')
         open (unit = 118, file = 
     &      'chains/dnde_best.dat',
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
         amptemp(1:sdim) = ampbest(1:sdim)
         alphap0 =  ampbest(sdim -1)
         normprim0 = ampbest(sdim)
         flagnorm =1
         call setprie_exp( alphap0,normprim0,flagnorm)
         call setsecep_exp   ! calculate the secondary electron and positron 
         call addprisece 
         write(*,*) 'chi2 = ', chi2Function(ampbest)

         call   priep_kra(alphap0,normprim0, flagnorm)

            
       
      end 




