module DataLoad
use params
implicit none
   public
   private i, j, alive,filename,fileid
    
   integer i, j
   logical alive
   character(len=79) :: filename
   character(len=79) :: datadir ='eptestF/data/'
   integer fileid
  
   integer, parameter:: ratelnum = 5929
   integer, parameter:: rateenglnum = 77
   real(kind=8):: PositronRate(ratelnum,3) , ElectronRate(ratelnum,3)
   real(kind=8):: RateEngList(rateenglnum)
   real(kind=8):: PositronRateMatrix(rateenglnum, rateenglnum)
   real(kind=8):: ElectronRateMatrix(rateenglnum, rateenglnum)
   integer, parameter:: AMSprlnum = 96
   real(kind=8):: AMSprotonE( AMSprlnum), AMSprotonF(AMSprlnum),     &
      AMSprotonErru(AMSprlnum), AMSprotonErrd(AMSprlnum)

   integer, parameter:: AMSeplnum = 59
   real( kind =8):: AMSepE(AMSeplnum) ,AMSepF(AMSeplnum),  &
      AMSepErru(AMSeplnum), AMSepErrd(AMSeplnum)
   real(kind=8):: AMSepEmin(AMSeplnum), AMSepEmax(AMSeplnum)

   integer, parameter:: AMSelnum = 63
   real( kind =8):: AMSelectronE(AMSelnum) ,AMSelectronF(AMSelnum),  &
      AMSelectronErru(AMSelnum), AMSelectronErrd(AMSelnum)
   real(kind=8):: AMSelectronEmin(AMSelnum),  &
      AMSelectronEmax(AMSelnum)

   integer, parameter:: AMSBClnum = 18
   real( kind =8):: AMSBCE(AMSBClnum),AMSBC(AMSBClnum),              &
      AMSBCErrd(AMSBClnum), AMSBCErru(AMSBClnum)
   integer, parameter:: AMSpolnum = 62 

   real( kind =8) :: AMSpositronE(AMSpolnum),AMSpositronF(AMSpolnum),&
      AMSpositronErru(AMSpolnum), AMSpositronErrd(AMSpolnum)
   real(kind=8):: AMSpositronEmin(AMSpolnum),                        &
      AMSpositronEmax(AMSpolnum)

   integer, parameter:: AMSperlnum = 65
   real(kind=8):: AMSperE( AMSperlnum), AMSperEmin(AMSperlnum),      &
      AMSperEmax(AMSperlnum), AMSper(AMSperlnum) ,                   &
      AMSperErru(AMSperlnum), AMSperErrd(AMSperlnum)
   
   
   real(kind=8):: y2AMSprotonF
   real(kind=8):: ElectronSource0(rateenglnum) 
   real(kind=8):: y2ElectronSource0(rateenglnum) 
   real(kind=8):: PositronSource0(rateenglnum)
   real(kind=8):: y2PositronSource0(rateenglnum)
   real(kind=8),parameter:: PI = 3.14159265359
   

   real(kind=8)::Sece_AMSelectron(AMSelnum)
   real(kind=8)::Sece_AMSpositron(AMSpolnum)
   real(kind=8)::Sece_AMSep(AMSeplnum)
   real(kind=8)::Secpo_AMSelectron(AMSelnum)
   real(kind=8)::Secpo_AMSpositron(AMSpolnum)
   real(kind=8)::Secpo_AMSep(AMSeplnum)


   real(kind=8)::AMSelectronSpcV(stepnum,AMSelnum)
   real(kind=8)::AMSpositronSpcV(stepnum,AMSpolnum)
   real(kind=8)::AMSepSpcV(stepnum,AMSeplnum)
   
contains

   subroutine SecSourceIni()
      call epRateLoad
      call AMSprotonLoad
      call SecSourceStep0
   end subroutine

   subroutine AMSDataIni()
      call AMSelectronLoad
      call AMSAllelectronLoad
      call AMSBtoCLoad
      call AMSPositronLoad
      call AMSperLoad
   end subroutine



   subroutine AMSperLoad()
! AMS positron to electron ratio
      real(kind=8):: temperr1, temperr2
      filename =trim(datadir)//'AMS02per2012.dat'
      fileid =180
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSperlnum
         read( fileid,*) AMSperE(i),AMSperEmin(i),AMSperEmax(i),     & 
            AMSper(i), temperr1,temperr2
         AMSperErru(i)= sqrt( temperr1**2 + temperr2**2)
         AMSperErrd(i) =AMsperErru(i)
      enddo
      close(fileid)
      
   end subroutine



   subroutine AMSelectronLoad()
! Unit : GeV^-1 s^-1 sr^-1 cm^-2   
      filename =trim(datadir)//'ElectronFlux_AMS02_fromPlots.dat'
      fileid =170
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSelnum
         read( fileid,*) AMSelectronE(i),AMSelectronF(i)
         AMSelectronF(i)= AMSelectronF(i)/(AMSelectronE(i)**3)*1.d-4
         if ( AMSelectronE(i) < 60.d0) then
            AMSelectronErrd(i)= 4.d-2 * AMSelectronF(i)
            AMSelectronErru(i)= 4.d-2 * AMSelectronF(i)
         else if ( AMSelectronE(i) < 1.0d2) then
            AMSelectronErrd(i)= 5.d-2 * AMSelectronF(i)
            AMSelectronErru(i)= 5.d-2 * AMSelectronF(i)
         else
            AMSelectronErrd(i)= 7.d-2 * AMSelectronF(i)
            AMSelectronErru(i)= 7.d-2 * AMSelectronF(i)
         endif
      enddo

      AMSelectronEmin(1)=sqrt(AMSelectronE(1)**3 / AMSelectronE(2)) 
      do i = 1, AMSelnum-1
         AMSelectronEmin(i+1)=sqrt(AMSelectronE(i)*AMSelectronE(i+1)) 
         AMSelectronEmax(i)= AMSelectronEmin(i+1)
      enddo
      AMSelectronEmax(AMSelnum)=sqrt(AMSelectronE(AMSelnum)**3 /   &
         AMSelectronE(AMSelnum-1) )
      !do i = 1, AMSelnum
         !write(*,*) AMSelectronEmin(i), AMSelectronE(i),             &
            !AMSelectronEmax(i)
      !enddo
      close(fileid)
      
   end subroutine

   

   subroutine AMSPositronLoad()
      filename =trim(datadir)//'PositronFlux_AMS02_fromPlots.dat'
      fileid =160
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSpolnum
         read( fileid,*) AMSpositronE(i),AMSpositronErrd(i),  &
            AMSpositronF(i),AMSpositronErru(i)
         AMSpositronErru(i)=abs(AMSpositronErru(i)-AMSpositronF(i))/ & 
            (AMSpositronE(i)**3) * 1.d-4
         AMSpositronErrd(i)=abs(AMSpositronErrd(i)-AMSpositronF(i))/ &
            (AMSpositronE(i)**3) * 1.d-4
         AMSpositronF(i)= AMSpositronF(i)/(AMSpositronE(i)**3)*1.d-4
      enddo

      AMSpositronEmin(1)=sqrt(AMSpositronE(1)**3 / AMSpositronE(2)) 
      do i = 1, AMSpolnum-1
         AMSpositronEmin(i+1)=sqrt(AMSpositronE(i)*AMSpositronE(i+1)) 
         AMSpositronEmax(i)= AMSpositronEmin(i+1)
      enddo
      AMSpositronEmax(AMSpolnum)=sqrt(AMSpositronE(AMSpolnum)**3 /   &
         AMSpositronE(AMSpolnum-1) )
      !do i = 1, AMSpolnum
         !write(*,*) AMSpositronEmin(i), AMSpositronE(i),             &
            !AMSpositronEmax(i)
      !enddo

      close(fileid)
      
   end subroutine


   subroutine AMSBtoCLoad()
      filename =trim(datadir)//'BoCRatio_AMS02_fromPlots.dat'
      fileid =150
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSBClnum
         read( fileid,*) AMSBCE(i),AMSBCErrd(i),AMSBC(i),AMSBCErru(i)
         AMSBCErru(i)=abs(AMSBCErru(i)-AMSBC(i))
         AMSBCErrd(i)=abs(AMSBCErrd(i)-AMSBC(i))
      enddo
      close(fileid)
      
   end subroutine



   subroutine AMSAllelectronLoad()
! Unit : GeV^-1 s^-1 sr^-1 cm^-2   
      filename =trim(datadir)//'AllElectronFlux_AMS02_fromPlots.dat'
      fileid =140
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSeplnum
         read( fileid,*) AMSepE(i),AMSepErrd(i),  &
            AMSepF(i),AMSepErru(i)
         AMSepErru(i)=abs(AMSepErru(i)-AMSepF(i))/ & 
            (AMSepE(i)**3) * 1.d-4
         AMSepErrd(i)=abs(AMSepErrd(i)-AMSepF(i))/ &
            (AMSepE(i)**3) * 1.d-4
         AMSepF(i)= AMSepF(i)/(AMSepE(i)**3)*1.d-4
      enddo
      AMSepEmin(1)=sqrt(AMSepE(1)**3 / AMSepE(2)) 
      do i = 1, AMSeplnum-1
         AMSepEmin(i+1)=sqrt(AMSepE(i)*AMSepE(i+1)) 
         AMSepEmax(i)= AMSepEmin(i+1)
      enddo
      AMSepEmax(AMSeplnum)=sqrt(AMSepE(AMSeplnum)**3 /   &
         AMSepE(AMSeplnum-1) )
      do i = 1, AMSeplnum
         write(*,*) AMSepEmin(i), AMSepE(i), AMSepEmax(i)
      enddo
      close(fileid)
      
   end subroutine

   




   subroutine epRateLoad ()
      integer:: fileidlist(2) = (/110, 120/)
      character(len=79) :: filenamelist(2)
      integer::cti

      filenamelist(1) = trim(datadir)// 'electrons_xsec.dat'
      filenamelist(2) = trim(datadir)// 'positrons_xsec.dat'
      do i = 1, 2
         inquire( file = filenamelist(i), exist=alive)
         if(.not. alive) then
            write(*,*) trim(filenamelist(i)),"doesn't exist."
            stop
         end if
         open ( fileidlist(i), file = filenamelist(i))
         read ( fileidlist(i) ,*) 
         do j = 1, ratelnum
            if ( i .eq. 1) then 
              read( fileidlist(i),*) ElectronRate(j,1), &
               ElectronRate(j,2), ElectronRate(j,3)
            else if ( i .eq. 2 )  then
              read( fileidlist(i),*) PositronRate(j,1), &
               PositronRate(j,2), PositronRate(j,3)
            endif
         enddo
         close( fileidlist(i))
      enddo
      do i = 1, rateenglnum
         RateEngList(i) = ElectronRate(i ,2)
      enddo
      do i = 1, rateenglnum
         do j = 1, rateenglnum
            cti = (i-1)*rateenglnum+ j
            
            PositronRatematrix(i,j)= PositronRate(cti,3)*1.d-24
            ! 1.d-24 transfrom barn to cm^2
            ElectronRatematrix(i,j)= ElectronRate(cti,3)*1.d-24
         enddo
      enddo
      
   end subroutine
   
   subroutine AMSprotonLoad()
      filename = trim(datadir)//'ProtonsFlux_AMS02_fromPlots.dat'
      fileid =130
      inquire ( file = filename, exist= alive)
      if(.not. alive) then
         write(*,*) trim(filename),"doesn't exist."
         stop
      end if
      open ( fileid, file = filename)
      read ( fileid ,*) 
      do i = 1, AMSprlnum
         read( fileid,*) AMSprotonE(i),AMSprotonErrd(i),  &
            AMSprotonF(i),AMSprotonErru(i)
         AMSprotonErru(i)= abs(AMSprotonErru(i) - AMSprotonF(i))
         AMSprotonErrd(i)=abs(-AMSprotonErrd(i) + AMSprotonF(i))
             
      enddo
      close(fileid)
      
      call spline( AMSprotonE, AMSprotonF,AMSprlnum,0.d0,0.d0, &
         y2AMSprotonF)

   end subroutine
   


   function AMSProton_Flux( eng)
! Unit GeV^-1 cm^-2 sr^-1 s^-1
      real(kind=8) :: AMSProton_Flux,eng,res
      
      AMSProton_Flux =0.
      if ( eng < AMSprotonE(1) ) then
         AMSProton_Flux = AMSprotonF(1)
         return
      elseif ( eng > AMSprotonE( AMSprlnum)) then
         AMSProton_Flux = 0.d0
         return
      end if 

      call splint(AMSprotonE,AMSprotonF,y2AMSprotonF,AMSprlnum,   &
        eng,res)
      AMSProton_Flux = res

   end function
   
   subroutine SecSourceStep0()
   ! calucate the coordinate indpendent part of source 
   ! Q / n(x)
   ! To calculate the Sec Source ( DO NOT FORGET to initilize
   ! the function 
      ElectronSource0 = 0.d0
      PositronSource0 = 0.d0
      do i = 1, rateenglnum
         do j = 1, rateenglnum
            ElectronSource0(i) = ElectronSource0(i) +            &  
               AMSProton_Flux( RateEngList(j))* 4. *PI *         &
               ElectronRateMatrix(i,j) * RateEngList(j)*0.2
            PositronSource0(i) = PositronSource0(i) +            &  
               AMSProton_Flux( RateEngList(j))* 4. *PI *         &
               ElectronRateMatrix(i,j) * RateEngList(j)*0.2 
         enddo
      enddo
      call spline(RateEngList, ElectronSource0,rateenglnum,0.d0, &
         0.d0, y2ElectronSource0)
      call spline(RateEngList, PositronSource0,rateenglnum,0.d0, &
         0.d0, y2PositronSource0)
   end subroutine
  
   function SecElectronSource(eng, rkpc, zkpc)
      use HydroDensity
   ! I think Proton flux is constant every where,
   ! how much we can trust the result
   ! Unit: GeV^-1 cm^-3 s^-1
      real(kind=8)::SecElectronSource 
      real(kind=8):: eng, rkpc, zkpc
      real(kind=8):: res, res2
      SecElectronSource =0.d0
      if ( eng .lt. RateEngList(1) .or. eng .gt.                 &
        RateEngList(rateenglnum) ) then
         write(*,*) 'Warnings!!!:SecElectron energy exceed the table'
         return
      endif
         
      call splint(RateEngList,ElectronSource0,y2ElectronSource0, &
        rateenglnum, eng,res)
      call nTotal_gal( rkpc,zkpc, res2)
      SecElectronSource = res * res2 
   end function

   function SecPositronSource(eng, rkpc, zkpc)
      use HydroDensity
   ! I think Proton flux is constant every where,
   ! how much we can trust the result
   ! Unit: GeV^-1 cm^-3 s^-1
      real(kind=8)::SecPositronSource
      real(kind=8):: eng, rkpc, zkpc
      real(kind=8):: res, res2
      SecPositronSource =0.d0
      if ( eng .lt. RateEngList(1) .or. eng .gt.                 &
        RateEngList(rateenglnum) ) then
         write(*,*) 'Warnings!!!:SecPositron energy exceed the table'
         return
      endif
         
      call splint(RateEngList,PositronSource0,y2PositronSource0, &
        rateenglnum, eng,res)
      call nTotal_gal( rkpc,zkpc, res2)
      SecPositronSource = res * res2 
   end function

end module
      
      
