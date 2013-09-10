module params
implicit none

! Model Parameters
        ! total dimension   spectrum aplitude+ 
        ! shift of Fermidata + index and  
        !   normalization of primary electron
        integer sdim
        parameter(sdim=10)

        !dimensionality of bin amplitude
        integer sdim1
        parameter(sdim1=7)

        !dimensionality of loglog amplitude = sdim1 -1
        integer sdimlog
        parameter(sdimlog=6)
         
        ! additional dimension 
        integer sdim2
        parameter(sdim2=3)

        ! the number of step function for tabulation,
        ! the number should  equal to sdim1
        integer stepnum 
        parameter(stepnum=7)
        double precision Famp(stepnum)

        ! chi2 energy cut
         double precision EngCut   
         parameter( EngCut = 10.d0)

        ! prefix for saved files
        character*100 prefix_f 
        parameter(prefix_f="chains/DmWoFermiLoglog7_2p7_0p3_")
        
        
        ! mass step for DM spectrum: Normalize at 1 TeV 
        ! Can be used for Pulsar also
        double precision msteptb
        parameter ( msteptb = 0.3d0)
        !parameter ( msteptb = 0.30d0)

        ! engergy step for pulsar spectrum 
        double precision engsteptb
        parameter ( engsteptb = 0.3d0)
        !parameter ( engsteptb = 0.30d0)
       
        !  highest log10 energy 
        double precision hlog10
        parameter(hlog10 = 2.7d0)  ! 4 corrsponds to 1.e4 GeV

 

      
        ! PAMELA electron data
        integer pamelanum
        parameter(pamelanum=39)
        double precision pamelaE(39), pamelaerr(39)
        double precision pamelaspcV(stepnum,39)
        double precision pamelaflux(39)

        ! Fermi data
        integer flnum
        parameter(flnum=42)
        double precision FEmid(42), Ferrp(42), Ferrm(42)
        double precision FspcV(stepnum,42)
        double precision Fflux(42)
        double precision Ffluxy2(42),Ferry2(42)

        ! AMS data
        integer amsnum
        parameter(amsnum=65)
        double precision AMSE(65), AMSerr(65)
        double precision AMSspcV(stepnum,65)
        double precision AMSpr(65)

        ! HESS data
        integer hessnum
        parameter(hessnum=9)
        double precision hessE(9), hesserr(9)
        double precision hessspcV(stepnum,9)
        double precision hessflux(9)
 
        !background elctron and positron at Fermi and AMS
        ! PAMELA and HESS
        double precision bg_electron(42), bg_positron(42) ! Fermi
        double precision bg_e_AMS(65), bg_p_AMS(65)!dm_AMS(65) ! AMS
        double precision bg_e_pamela(39) ! pamela electron
        double precision bg_e_hess(9), bg_p_hess(9)! hess 
        double precision bg_prie_F(42),bg_sece_F(42),bg_secp_F(42)
        double precision bg_prie_A(65),bg_sece_A(65),bg_secp_A(65)
        double precision bg_prie_P(39),bg_sece_P(39),bg_secp_P(39)
        double precision bg_prie_H(9) ,bg_sece_H(9) ,bg_secp_H(9)
        
        !tablulate primary electron index from 2. to 4.
        integer stepalpha 
        parameter(stepalpha=21)
        double precision prim_alpha(stepalpha)
        double precision prim_F_table(stepalpha,42) 
        double precision prim_A_table(stepalpha,65) 
        double precision prim_P_table(stepalpha,39) 
        double precision prim_H_table(stepalpha,9) 
        double precision prim_F_y2(stepalpha,42) 
        double precision prim_A_y2(stepalpha,65) 
        double precision prim_P_y2(stepalpha,39) 
        double precision prim_H_y2(stepalpha,9) 

        !tablulate dark matter spectrum index 
        integer stepalphalog 
        parameter(stepalphalog=61)
         ! alphaloglist set in main file
        double precision alphaloglist(stepalphalog)
        double precision dm_F_table(sdimlog,42,stepalphalog) 
        double precision dm_A_table(sdimlog,65,stepalphalog) 
        double precision dm_P_table(sdimlog,39,stepalphalog) 
        double precision dm_H_table(sdimlog,9,stepalphalog) 
        double precision dm_F_y2(sdimlog,42,stepalphalog) 
        double precision dm_A_y2(sdimlog,65,stepalphalog) 
        double precision dm_P_y2(sdimlog,39,stepalphalog) 
        double precision dm_H_y2(sdimlog,9,stepalphalog) 

        ! amplitude  Cube
         double precision amptemp(sdim)
           
             
      	!no. of modes to generate
      	!integer sModes 
	!parameter(sModes=2)
      
      	!width of the Gaussian profile of each ring
	!double precision sw(sModes)
	!data sw /0.1d0,0.1d0/
      
      	!width of the rings
      	!double precision sr(sModes)
	!data sr /2.d0,2.d0/
      
      	!Center of the rings. 
      	!Centers in the 1st dimensions are set here while in 
      	!all the other dimensions, they are set to 0 in main.f90
      	!double precision sc(sModes,sdim)
	!data sc(1,1),sc(2,1) /-3.5d0,3.5d0/
      
      	!priors on the parameters
      	double precision spriorran(sdim,2)
     
        !sigma of the Gaussian 
!        double precision sigma(sdim)
 


! Parameters for Nested Sampler
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.true.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.false.)
	
      	!max no. of live points
      	integer nest_nlive
	parameter(nest_nlive=10000)
      
      	!tot no. of parameters, should be sdim in most cases but if you need to
      	!store some additional parameters with the actual parameters then
      	!you need to pass them through the likelihood routine
	integer nest_nPar 
	parameter(nest_nPar=sdim)
      
      	!seed for nested sampler, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-5)
      
      	!evidence tolerance factor
      	double precision nest_tol 
      	parameter(nest_tol=0.5)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
      	parameter(nest_efr=0.5d0)
      
      	!root for saving posterior files
      	character*100 nest_root
	parameter(nest_root='chains/epchi2-')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=1000)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=20)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	parameter(nest_nClsPar=3)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	parameter(nest_resume=.false.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.true.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1d0))
      
      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	!parameters to wrap around (0 is F & non-zero T)
	integer nest_pWrap(sdim)
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================


end module params
