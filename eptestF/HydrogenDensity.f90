module HydroDensity
implicit none
   private i, j
   integer i, j
   
contains
   subroutine nH2_gal( r, z, res)  
   !res-> density of molecular hydrogen
      real( kind =8):: rlist(18) = (/0.00, 2.25, 2.75, 3.25,3.75,&
         4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25,   &
         8.75, 9.25, 9.75,10.25 /)
      double precision:: y(18) = (/ 0.0, 1.5, 3.3, 5.8, 5.5, 8.4,  &
         9., 9.6, 8.6, 9.1, 7.9, 9.2, 7.7, 5., 3.6, 4.8, 1.7, 0./) 
      double precision:: z0(18) = (/ 0.039,0.039,0.036,0.000,-.008,&
         0.001,-.010,-.001,-.004, -.019,-.022,-.014,-.009,-.004, &
         0.013,-.004,-.020,-.020 /)
      double precision:: zh(18) = (/ 0.077,0.077,0.080,0.061,0.065,&
         0.071,0.072,0.082,0.083, 0.073,0.063,0.058,0.072,0.080, &
         0.066,0.023,0.147,0.147 /)
      double precision r, z, res
      double precision H2toCO
      double precision fr, fz0, fzh 
      double precision kpc
      kpc = 3.08568d21
      
      H2toCO = 1.9*1.e20
      if ( r .gt. rlist(18) )  then
         res = 0.
         return
      end if
      i =1
      do while (  rlist(i+1) <= r  )
         i=  i+1
      end do
      
      fr  = y(i) + ( y(i+1) - y(i) ) / ( rlist(i+1)-rlist(i)) *  & 
         (r- rlist(i))  
      fz0 = z0(i) + ( z0(i+1) - z0(i))/ ( rlist(i+1)-rlist(i))*  &
         ( r - rlist(i)) 
      fzh = zh(i) + ( zh(i+1) - zh(i))/ ( rlist(i+1)-rlist(i))*  &
         ( r - rlist(i)) 
      res = fr * exp( - log(2.)* ( (z - fz0)/fzh)**2 ) * H2toCO
      res = res/ kpc
   end subroutine

   subroutine nHI_gal ( rkpc, zkpc, res)
      double precision rkpc, zkpc, res
      double precision:: r(30)= (/0., 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,&
         5.0, 5.5, 6., 6.5, 7.0, 7.5, 8.0, 8.5, 9., 9.5,10.,10.5,&
         11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0 /)
      double precision:: y(30)= (/.10, .13, .14, .16, .19, .25,.30,&
         .33, .32, .31, .30, .37, .38, .36, .32, .29, .38, .40,  & 
         .25,.23,.32, .36, .32, .25, .16, .1, .09, .08, .06,.00 /)
      double precision:: fr, fz, fz1 =0., fz2 =0., r1, r2 
      double precision y1, y2 
      double precision:: nGB =0.33, nDL  = 0.57
      double precision:: A1 = 0.395, z1 = 0.212/2.
      double precision:: A2 = 0.107, z2 = 0.530/2.
      double precision:: B = 0.064, zh = 0.403
      r2 = r(30)
      y2 = y(30)
      i =1
      do while (  r(i+1) <= rkpc  )
         i=  i+1
      end do
      
      r1  = ( r(i) + r(i+1) )/2. 
      y1= y(i)
      if ( rkpc < r1) then
         if ( i>1) then
            r2 = ( r(i-1) +r(i) )/2. 
            y2 = y(i-1)
         else
            r2 = r(1)
            y2 = y(1) 
         end if
      else if ( i < 29) then
         r2 = (r(i+1) + r(i+2))/2.
         y2 = y(i+1)
      end if

   ! interpolation in r
      fr = y1 +( y2-y1) / ( r2-r1)*( rkpc - r1)
   
      r2 = (r(29) + r(30))/2.
   ! extrapolation in r
      if ( rkpc > r2) fr = y(29)* exp( - ( rkpc - r2)/3.)
      
   ! Calucation of z-dependence
      if ( rkpc < 10.) then
         fz1 = A1 * exp( -log(2.)* (zkpc /z1)**2) + A2* exp(     &
            -log(2.)*(zkpc/z2)**2) + B*exp(-abs(zkpc)/zh)
      end if
      if (rkpc >8.) then
         fz2 = nDL*exp( - ( zkpc/ ( 0.0523*exp(0.11*rkpc)))**2) 
      end if
      if (rkpc <= 8.) then
         fz = fz1
      else 
         if ( rkpc >= 10.) then
            fz= fz2
         else 
            fz = fz1 + ( fz2- fz1)/ 2.* ( rkpc-8.)
         end if
      end if
      
      res = fz*fr/nGB
      
   end subroutine       



   subroutine nHII_gal( rkpc, zkpc, res)
      double precision rkpc, zkpc, res
      double precision:: fne1 = 0.025, H1 = 1.0 , A1 = 20.0
      double precision:: fne2 = 0.200, H2 =0.15 , A2 = 2.
      double precision:: r2 = 4.0
      
      double precision ne1, ne2
      ne1 = fne1 * exp( -abs(zkpc)/H1)* exp(-(rkpc/A1)**2)
      ne2 =fne2 * exp(-abs(zkpc)/H2)* exp(-((rkpc-r2)/A2)**2)
      res = ne1 + ne2
   end subroutine

   subroutine nTotal_gal( rkpc, zkpc, res)
      double precision rkpc, zkpc, res
      double precision n1, n2, n3
      call nH2_gal( rkpc, zkpc, n1)
      !write(*,*) 'n1',n1
      call nHI_gal( rkpc, zkpc, n2)
      !write(*,*) 'n2',n2
      call nHII_gal( rkpc, zkpc, n3)
      !write(*,*) 'n3',n3
      res = n1 * 2. + n2 + n3
   end subroutine

end module
      
      
