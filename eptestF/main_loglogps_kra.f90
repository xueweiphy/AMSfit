      include 'user-dummy-eptest.h'
      include 'user-not-dummy-eptest.h'

      program main
        use nestwrapper
        use params






ccc  Use Pulsar to explain data
ccc  I am not sure how to switch between pulsar and secondary electron 
ccc  and positron calculation



ccc routine to find the best fit to Fermi, HESS, AMS, PAMELA
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


ccccc   from pulsar test file
        real*8 dsepdphidp_setpul,deltat
        character*100 popfile
      common/pulsarfluxc/ popfile,deltat


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
      integer flagoutr


      integer pulsardmflag
      common /compsdm/ pulsardmflag

      real*8 dsepfluxdmps

      call dsinit


      call readdatafiles

      binnumle  = stepnum
      pulsardmflag  = 1  ! 1 for pulsar , 0 for dark matter 

      

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
      bfieldmean=7.d0
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
c      k0halo=k0halo*(3.0856d21)**2/(1.d6*31556926d0)*1.d-27 ! 10^27 cm^2 s^-1
      k0halo =   29. !  10^27 cm^2 s-1  
                     !  this value provided by Daniele
      k0gasdisc=k0halo
      write(*,*) 'K0 : ',k0halo
      kdiffdeltalow=0.d0  ! spectral index below the break - dummy variable
      kdiffdelta=0.5d0    ! spectral index above the break
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


ccc        Use PAMELA data at 30.6 GeV to normalize the flux


      alphap0= 2.7d0       
      normprim0 =1.d0
      flagnorm =1  ! 1 is to reset the normalization
      do i = 1, stepalpha
         alphap0 = 1.9 + 0.1*real(i)
         prim_alpha(i) = alphap0
         call setprie_exp( alphap0,normprim0,flagnorm)
         prim_F_table(i,:) =  bg_prie_F
         prim_A_table(i,:) =  bg_prie_A
         prim_P_table(i,:) = bg_prie_P
         prim_H_table(i,:) =  bg_prie_H
      enddo

      call spline_prim
      alphap0 = 2.70
      call splint_prim(alphap0,normprim0)

      call setsecep_exp   ! calculate the secondary electron and positr    on 
      call addprisece


c      alphap0 =  2.4149952567455548
c      normprim0=  9.8754429395077459
c      call   priep_kra(alphap0,normprim0, flagnorm)




ccc   SHOULD I SWITCH TO spherical distribution
ccc      npaxidefault=.false.  ! link to spherical distribution ???


cccccccccccccccccccccccccccccccccccccc
ccc             Pulsar             ccc
cccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccc
ccc       change modeep 
ccccccccc   6--- step function
cccccccc    7--- theta function
cccccccc    8--- loglgo interpolation function
cccccccc    12--- for tabulate loglog interpolation

      modeep=6 ! Here you should NOT set modeep to 8, 
               ! since we use 6 or 7 to set scan region

      phiin=0.0d0   ! forced field solar modulation parameter 
      deltat=0.d0

ccc
ccc this is the catalogue for our analysis !!!!!!!!!!!!!!!!!
ccc TAKE IT FROM HERE
ccc
      popfile='eptestF/pulsar-cat130508.txt'


      do k =1, stepnum
         engvec(k) = 10.d0**( hlog10- engsteptb*(k-1))
         engvecmin(k) = 10.d0**( hlog10- engsteptb*(k))
         amplitudeR(k) = 0.d0
         logenglist(stepnum+1-k)= hlog10- engsteptb*(k-1)
      enddo

      
      feps =1.d-3
      prec = 1.d-2
c      call spectrumkra
c          call writekra

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc tabulate spectrum . mass from ~1 TeV to ~10 GeV
ccc   energy are at Fermi, AMS, PAMELA, HESS data points
ccc
       intflag = 0 
       call spectrumdata( intflag) 
    

       amptemp(1:stepnum) = 1.d1
       WRITE(*,*)  'mode 7 -----------------'  
            modeep=7
         call spectrumdata(0)
       DO k = 1  ,1  
         do i = 1, 42  
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
             k=1
      write(*,*) 'mode =8'
         do i = 1, 42  
    
       WRITE(*,*) FEmid(i), FspcV(1,i)  
         enddo

            modeep=7
        call spectrumdata( intflag) 
ccccccccccccc  using Fermi Set MAX amplitude


      do k = 1, sdim1
         amptemp1=10.d10
         amptemp2=10.d20
         i=7
         do j = 6 , 42
            fluxtemp= Fflux(j)+ 4.*Ferrp(j) -
     &        bg_electron(j)-bg_positron(j)
c         write(*,*) k, j , 'flux residual', fluxtemp
            if ( fluxtemp > 0 .and. Fspcv(k,j) > 0.)
     &        amptemp2 =   fluxtemp /(2.* FspcV(k,j)  )
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
          if ( fluxtemp>0 .and. AMSspcv(k,j) > 0.)
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





ccc Begin TEST loglog interpolation Table 

      normlist_test= stepampmax  
      write(*,*) 'normlist_test'
      write(*,*) normlist_test
      modeep = 8 
      normlist_test =1.
      amptemp(1:sdim1) = normlist_test
      
c      call spectrumdm_data(0,0)  !spectrumdm also can be used for pulsar
                                 ! for modeep =10
      call spectrumdata(0)
      flux_F_test2 = FspcV(1,:)
      write(*,*) 'Flux_F_test2'
      write(*,*) flux_F_test2
      write(*,*) amplitudeR

      modeep = 12
      mstep1 = msteptb  
      do i = 1, 61
         alphaloglog  = (real(i)-31.0)*0.1 /mstep1
         alphaloglist(i) = alphaloglog
         listalphacom = alphaloglog
         call spectrumdm_data(1,0)
         dm_A_table(:,:,i)=AMSspcV(1:sdimlog,:)
         dm_F_table(:,:,i)=FspcV(1:sdimlog,:)
         dm_P_table(:,:,i)=pamelaspcV(1:sdimlog,:)
         dm_H_table(:,:,i)=hessspcV(1:sdimlog,:)
      enddo

      call spline_spcloglog
      call alpha_loglog(normlist_test,alphalog0,flagoutr)
      write(*,*) 'alphalog' 
      write(*,*) alphalog0 
      write(*,*) 'normlist' 
      write(*,*) normlist_test 
      flux_F_test =0.d0
      do i = 1 , sdimlog
        flux_F_test = flux_F_test + normlist_test(i)*FspcV(i,:)
      enddo
      open (unit = 102, file = 'Root/spectrum_F_test.dat',
     &  form='formatted' , status='Replace')

      write(*,*) ' write spectrum to spectrum_stepf'
      write(102,1020) FEmid
      write(102,1020) flux_F_test2
      write(102,1020) flux_F_test

      listalphacom(:sdimlog) = alphalog0
      listalphacom(sdimlog+1) = 0.d0
      call spectrumdm_data(0,0)
      flux_F_test =0.d0
      do i = 1 , sdimlog
        flux_F_test = flux_F_test + normlist_test(i)*FspcV(i,:)
      enddo
      write(102,1020) flux_F_test

1020  FORMAT(42(1x,e14.8) )

      close(102)




ccc TEST end




         !setting priors
   
        do i = 1, sdim1
           spriorran(i,1)= 0.
           spriorran(i,2)= stepampmax(i) *1e3
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
         
         
         modeep =12

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
         amptemp(1:sdim) = ampbest(1:sdim)
         alphap0 =  ampbest(sdim -1)
         normprim0 = ampbest(sdim)
         flagnorm =1
c         call setprie_exp( alphap0,normprim0,flagnorm)
         call splint_prim( alphap0,normprim0)
         call setsecep_exp   ! calculate the secondary electron and positron 
         call addprisece 
         intflag =0
         if ( modeep .eq. 10) then
            call spectrumdm_data(intflag,1)
         else if ( (modeep .eq. 11) .or. (modeep .eq. 12)) then
            call alpha_loglog(ampbest(1:sdim1),alphalog0,flagoutr)
            call splint_spcloglog(ampbest(1:sdim1),alphalog0)
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
            if (ppin .lt. hamwimp) then
               spcvec(k,i) =   dsepfluxdmps(ppin)
            else 
               spcvec(k,i) =   0.
            endif
         enddo
      enddo



      call writekra
      call writedata
      call   priep_kra(alphap0,normprim0, flagnorm)




            
       
       end 




