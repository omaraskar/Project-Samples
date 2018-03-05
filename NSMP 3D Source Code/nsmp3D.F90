!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	    MAIN PROGRAM                              !
!       Solution of the Navier-Stokes Multi-phase equations           !
!                      Miguel Angel Uh Zapata                         !
!                  Last modification: Jul 2014                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!
      PROGRAM nsmp3D
!
!---------------------------------------------------------------------!   
!                                                                     !
!     This program find the solution of the Navier-Stokes equations   !
!     corresponding to the free surface two-phase flow problem in 3D. !
!     In this program all the defined and common variables are modi-  !
!     fied.                                                           ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  | nstep       | Current number of time step                     |  !
!  | tcpu        | cpu time of the simulation                      |  !  
!  | tt(1:2)     | user and system simulation time                 |  !
!  | nt          | Approx number of time steps to display time     |  !
!  | time0       | Initial time of the current simulation          |  !
!  | usepas      | Elapsed time for one point & one time step      |  !  
!  | idhr        | Number of hours of the simulation               |  !
!  | idmin       | Number of minutes of the simulation             |  !  
!  | idsec       | Number of seconds of the simulation             |  !                
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |  STRUCTURAL:                                                  |  !  
!  |                      ° common.mpf                             |  !
!  |                      ° cppdefs.h                              |  !
!  |                      ° definition_variables                   |  !
!  |                      ° alloc_variables                        |  !
!  |                      ° dealloc_variables                      |  !
!  |                      ° interfaces.F90                         |  !
!  |                                                               |  ! 
!  |  INITIALIZATION:                                              |  ! 
!  |                      ° input.F90                              |  !
!  |                      ° input_parameters.F90                   |  !
!  |                      ° geometry.F90                           |  !
!  |                      ° initial.F90                            |  !
!  |                                                               |  !
!  |  TIME LOOP UPDATES:                                           |  ! 
!  |                      ° LoopTimeInitial.F90                    |  !
!  |                      ° hydro.F90                              |  !
!  |                      ° LoopTimeUpdate.F90                     |  !
!  |                                                               |  !
!  |  RE-START & SAVING:                                           |  !
!  |                      ° outsavtec.F90                          |  !
!  |                      ° restart_out.F90                        |  !  
!  |                      ° restart_in.F90                         |  !     
!  |                                                               |  !
!  |_______________________________________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                            Definitions                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |            Keys, subroutines and parameters            |
!     |________________________________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE variables 
            USE geometry
            USE interfaces
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE variables 
            USE geometry 
            USE interfaces
            implicit none
#     endif
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|

      real,dimension(2) :: tt
      real ::tcpu,MAXtcpu
      real*8 :: nt,usepas
      real*8 :: time0
      integer:: idmin,idhr,idsec
!     --------------------------
      integer:: nstep  
      integer:: nRK
!     --------------------------
      character*80 title

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                            '      
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'------------------------------------------------------------' 
      print*,'            __________________________________              '
      print*,'           |                                  |             '
      print*,'           |    ***  PROGRAM: NSMP3D   ***    |             '
      print*,'           |__________________________________|             '
      print*,'                                                            '
      print*,'------------------------------------------------------------'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
      print*,'    ___________________________________________________     '
      print*,'   |                                                   |    '
      print*,'   |                  INITIALIZATION                   |    '
      print*,'   |___________________________________________________|    '
      print*,'                                                            '
#     endif
!     =============== END ================    
!     ====================================

!*********************************************************************!
!                                                                     !
!                         0) Parallelization                          !
!                                                                     !
!*********************************************************************!

!      _____________________________________________________________ 
!     |      |                                                      |
!     | 0.0  |     Read: Number of cell-center and vertex points    |
!     |______|______________________________________________________|

      open(25,file='data.txt',status='old')
      read(25,*) title
      read(25,*) N_VERTglobal
      read(25,*) N_CELL0global
      close(25)
!     _____________________________________________________________ 
!    |      |                                                      |
!    | 0.1  |               Defition of parameters                 |
!    |______|______________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         NZ      = NZglobal
         N_VERT  = N_VERTglobal
         N_CELL0 = N_CELL0global
         N_CELL  = N_CELL0 + N_CELLghostMAX

         print*,'                                                       '
         print*,'       DOMAIN:                                         ' 
         print*,'          Number of elements  =',N_CELL0Global 
         print*,'          Number of vertices  =',N_VERTGlobal
         print*,'          Number of NZ points =',NZglobal-1
         print*,'                                                       '
#     endif
!     =============== END ================    
!     ====================================

!     _____________________________________________________________ 
!    |      |                                                      |
!    | 0.2  |                Parallel structure                    |
!    |______|______________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
!       __________________________________ 
!       Distribution         
        call initialisation_mpi
        call parallel_input_global
        call parallel_distribution
!       __________________________________ 
!       Topology
        call parallel_topology
        call parallel_neighbors
!       __________________________________  
!       Defition of parameters               
        call parallel_parameters
!       __________________________________ 
!       Parallel utilities 
        call parallel_index
        call parallel_shareindex
        call parallel_type      
#     endif
!     =============== END ================    
!     ====================================
     
!*********************************************************************!
!                                                                     !
!          1)           Variables & input data files                  !
!                                                                     !
!*********************************************************************!

!      _____________________________________________________________ 
!     |      |                                                      |
!     | 1.1  |                  Allocate variables                  |
!     |______|______________________________________________________|

      call alloc_variables
      call alloc_geometry
!      _____________________________________________________________ 
!     |      |                                                      |
!     | 1.2  |                    Input data                        |
!     |______|______________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
        call input(No_vp,No_cp,No_wb,No_hb,No_qb,No_sp,&
                   nbe,xv,yv,zbv) 
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel    
        call parallel_input_local(No_vp,No_cp,No_wb,No_qb,No_hb,No_sp,&
                                  nbe,xv,yv,zbv)
#     endif
!     =============== END ================    
!     ====================================
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.3  |                 Input parameters                       |
!     |______|________________________________________________________|

      call input_parameters(nbev,hv,h,                     &
                            No_vp,No_cp,No_wb,No_hb,No_qb, &
                            nbe,xv,yv,zbv)
!      _____________________________________________________________ 
!     |          |                                                  |
!     |    MPI   |              Test MPI Examples                   |
!     |__________|__________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel    
        !call TestMPInsmp2D(No_cp,nbe)
        !call TestMPInsmp3D(No_cp,nbe)
#     endif
!     =============== END ================    
!     ====================================


!*********************************************************************!
!                                                                     !
!                           2) Initialization                         !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 2.1  |            Geometry variables of the mesh              |
!     |______|________________________________________________________|
        
      call calcul_geometry(xc,yc,sig,dsig,No_cp,nbe,h, &
                           xv,yv,sigv,dsigv,No_vp,nbev,&
                           ic1tec,ic2tec,ic3tec)
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 2.2  |                Open time error file                    |
!     |______|________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         if (RunTimeExample.eq.1) then
            open(8100,file="../output/ErrorTime.dat",status='unknown')
         endif 
         if (RunFreeSurface.eq.1) then
            open(9100,file="../output/FS/FS_ErrTime.dat",status='unknown')
            open(9101,file="../output/FS/FS_SolTime1.dat",status='unknown')
            open(9102,file="../output/FS/FS_SolTime2.dat",status='unknown')
         endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         if (RunTimeExample.eq.1) then
            open(8100,file="../output/ErrorTime.dat",status='unknown')
         endif 
         if (RunFreeSurface.eq.1) then
            open(9100,file="../output/FS/FS_ErrTime.dat",status='unknown')
            open(9101,file="../output/FS/FS_SolTime1.dat",status='unknown')
            open(9102,file="../output/FS/FS_SolTime2.dat",status='unknown')
         endif
         ENDIF
#     endif
!     =============== END ================    
!     ====================================
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 2.3  |                 Initial conditions                     |
!     |______|________________________________________________________|

!     ___________________________________________
!     Initial time and step 

      time  = 0.0d0
      nstep = 0
!     ___________________________________________
!     Initial values of the main variables 

      call initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,   &
                   alphasn,usn,vsn,wsn,psn,viscos,rhos,   &
                   alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv, &
                   alphasv,usv,vsv,wsv,psv,viscosv,rhosv, &
                   etan,etav,Hpr,Hprv,                    &
                   Heaviside,mask,                        &
                   xct,yct,zct,                           &
                   xvt,yvt,zvt,                           &
                   xc,yc,sig,dsig,No_cp,nbe,              &
                   xv,yv,sigv,dsigv,No_vp,nbev,           &
                   h,hv)
!     ___________________________________________
!     Initial values of time test problems 

      if (RunTimeExample.eq.1) then
         call testTimeInitial(ufn,vfn,wfn,pfn,           &
                              ufv,vfv,wfv,pfv,           &
                              xc,yc,sig,dsig,No_cp,nbe,  &
                              xv,yv,sigv,dsigv,No_vp,nbev) 
      endif
!     ___________________________________________
!     Initial values of free surface problems 
      if (RunFreeSurface.eq.1) then
         call FS_testTimeInitial(ufn,vfn,wfn,pfn,             &
                                 ufv,vfv,wfv,pfv,             &
                                 xc,yc,sig,dsig,No_cp,nbe,    &
                                 xv,yv,sigv,dsigv,No_vp,nbev, & 
                                 Hpr,h,etan,                  &
                                 Hprv,hv,etav,                &
                                 xct,yct,zct,                 & 
                                 xvt,yvt,zvt,                 &
                                 No_sp) 
      endif

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 2.4  |               Re-start solution (n)                    |
!     |______|________________________________________________________|

!     ________________________________________________________
!     Save counter
      SaveCounter = 0

!     ________________________________________________________
!     Re-start time 

      if (IrestartIN.eq.1) then
!        ---------------------------------------------------
!        Read main variables      
         call restart_in(alphafn,ufn,vfn,wfn,pfn,  & 
                         alphasn,usn,vsn,wsn,psn,  &
                         zct)
!        ---------------------------------------------------
!        Update the free surface  
         do i=1,N_CELL
            etan(i)= zct(i,NZ)
            Hpr(i) = etan(i)+h(i)
         enddo        
!        ---------------------------------------------------
!        Display re-start initial time
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
         print*,'    ==================================================  '
         print*,'                  Re-start simulation                   '
         print*,'    ==================================================  '
         print*,'                                                        '       
         write(*,'(t17,a13,f10.3)') '       time = ', time
         write(*,'(t17,a13,f10.3)') '         dt = ', dt
         write(*,'(t17,a13,i8)')    'SaveCounter = ', SaveCounter
         print*,'                                                        '
         print*,'    ==================================================  '
         print*,'                                                        '
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
         IF (rang_topo.eq.0) THEN
         print*,'    ==================================================  '
         print*,'                  Re-start simulation                   '
         print*,'            WARNING!! Parallel not done yet             '
         print*,'    ==================================================  '
         print*,'                                                        '       
         write(*,'(t17,a13,f10.3)') '       time = ', time
         write(*,'(t17,a13,f10.3)') '         dt = ', dt
         write(*,'(t17,a13,i8)')    'SaveCounter = ', SaveCounter
         print*,'                                                        '
         print*,'    ==================================================  '
         print*,'                                                        '
         ENDIF
#        endif
!        =============== END ================    
!        ====================================
      endif

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 2.5  |             Save initial condition values              |
!     |______|________________________________________________________|

!     -------------------------------------------------------
!     Last time when the results were saved
      LastTimeSave = time

      if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
!        ---------------------------------------------------
!        Save Tecplot file cell centers       
         call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,              &
                            alphasv,usv,vsv,wsv,psv,rhosv,        &
                            xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save Tecplot file cell centers 
         call SavetecCenter(alphafn,ufn,vfn,wfn,pfn,              &
                            alphasn,usn,vsn,wsn,psn,rhos,         &
                            xct,yct,zct,No_cp,                    &
                            ic1tec,ic2tec,ic3tec)
      endif

!*********************************************************************!
!                                                                     !
!              3) Main structure: time step simulations               !
!                                                                     !
!*********************************************************************!

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         print*,'                                                        '
         print*,'    ___________________________________________________ '
         print*,'   |                                                   |'
         print*,'   |                TIME SIMULATIONS                   |'
         print*,'   |___________________________________________________|'
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         print*,'                                                        '
         print*,'   ==================================================== '
         print*,'                    MPI: TIME SIMULATIONS               '
         print*,'   ==================================================== '
         print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================    
!     ====================================

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.1  |          Time loop:  Initial values: (n+1)=(n)         |
!     |______|________________________________________________________|

      call LoopTimeInitial(etanp,                              &
                           alphafnp,ufnp,vfnp,wfnp,pfnp,       &
                           alphasnp,usnp,vsnp,wsnp,psnp,       &                                 
                           xct,yct,zct,                        & 
                           xvt,yvt,zvt,                        &    
                           etan,                               &  
                           alphafn,ufn,vfn,wfn,pfn,            &
                           alphasn,usn,vsn,wsn,psn,            &
                           xc,yc,sig,Hpr,h,                    &
                           xv,yv,sigv,Hprv,hv,                 &
                           No_cp,nbe)  

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.2  |         Time loop:  Beggining of the new time step     |
!     |______|________________________________________________________|
        
10    continue

      time  = time + dt
      nstep = nstep + 1

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifndef KeyDisplay
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(*,'(a8,i5)')  'Step : ',nstep
         write(*,'(a8,f8.4)')'Time = ',time
         print*,'  '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         write(*,'(a8,i5)')  'Step : ',nstep
         write(*,'(a8,f8.4)')'Time = ',time
         print*,'  '
         ENDIF
#     endif
!     =============== END ================    
!     ====================================
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________ 
!     |      |                                                 |
!     | 3.2.1|  Runge-Kutta: Initial known values: (known)=(n) |
!     |______|_________________________________________________|

      do i=1,N_CELL
!        ----------------------------
!        Free surface  
         eta(i)=etan(i)
         do k=1,NZ
!           -------------------------
!           Volume Fraction  
            alphaf(i,k) = alphafn(i,k)
            alphas(i,k) = alphasn(i,k)
!           -------------------------
!           Velocity
            uf(i,k) = ufn(i,k)
            vf(i,k) = vfn(i,k)
            wf(i,k) = wfn(i,k)
            us(i,k) = usn(i,k)
            vs(i,k) = vsn(i,k)
            ws(i,k) = wsn(i,k)
!           -------------------------
!           Pressure
            pf(i,k) = pfn(i,k)
            ps(i,k) = psn(i,k)
        enddo  
      enddo 

!      ________________________________________________________ 
!     |      |                                                 |
!     | 3.2.2|  Runge-Kutta: The two step loop                 |
!     |______|_________________________________________________|
      
      DO nRK=1,FinRK

!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        ifdef KeyDisplay
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
            print*,' '
            write(*,'(t15,a8,i3)')'RK step:',nRK  
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
            IF (rang_topo.eq.0) THEN
            print*,' '
            write(*,'(t15,a8,i3)')'RK step:',nRK  
            ENDIF
#        endif
!        =============== END ================    
!        ====================================
#        endif
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!         _______________________________________
!        |                                       |
!        |            Update variables           |
!        |_______________________________________|

!        ________________________________________
!        Calculate test examples (one step)  
         if (RunTestExample.eq.1) then	
            call testExamples(                                &
!                   -------------------------------------------
                    ufnp,ufn,uf,ufv,                          &
                    vfnp,vfn,vf,vfv,                          &
                    wfnp,wfn,wf,wfv,                          &
                    pfnp,pfn,pf,pfv,                          &
!                   -------------------------------------------
                    usnp,usv,                                 &
                    vsnp,vsv,                                 &
                    wsnp,wsv,                                 &
                    psnp,psv,                                 &
!                   -------------------------------------------
                    xc,yc,sig,dsig,No_cp,nbe,                 &
!                   --------------------------------------------
                    xv,yv,sigv,dsigv,No_vp,nbev,              &
!                   -------------------------------------------- 
                    Hpr,h,etan,                               &
!                   -------------------------------------------- 
                    Hprv,hv,etav,                             &
!                   --------------------------------------------
                    nRK)
         endif
!        ________________________________________
!        Calculate test time examples
         if (RunTimeExample.eq.1) then	
            call testTimeExamples(                            &
!                   -------------------------------------------
                    ufnp,ufn,uf,ufv,                          &
                    vfnp,vfn,vf,vfv,                          &
                    wfnp,wfn,wf,wfv,                          &
                    pfnp,pfn,pf,pfv,                          &
!                   -------------------------------------------
                    usnp,usv,                                 &
                    vsnp,vsv,                                 &
                    wsnp,wsv,                                 &
                    psnp,psv,                                 &
!                   -------------------------------------------
                    xc,yc,sig,dsig,No_cp,nbe,                 &
!                   --------------------------------------------
                    xv,yv,sigv,dsigv,No_vp,nbev,              &
!                   -------------------------------------------- 
                    Hpr,h,etan,                               &
!                   -------------------------------------------- 
                    Hprv,hv,etav,                             &
!                   --------------------------------------------
                    nRK)
         endif
!        ________________________________________
!        Free surface examples         
         if (RunFreeSurface.eq.1) then	
            if (RunTimeExample.eq.1) then
                print*, ' Error!!! RunTimeExample should be = 0'
                stop
            endif
            call FS_testTimeEqn(                              &
!                   -------------------------------------------
                    ufnp,ufn,uf,ufv,                          &
                    vfnp,vfn,vf,vfv,                          &
                    wfnp,wfn,wf,wfv,                          &
                    pfnp,pfn,pf,pfv,                          &
!                   -------------------------------------------
                    rhof,viscof,                              &
                    rhofv,viscofv,                            &
!                   -------------------------------------------
                    xc,yc,sig,dsig,No_cp,nbe,                 &
                    xv,yv,sigv,dsigv,No_vp,nbev,              &
!                   -------------------------------------------- 
                    Hpr,h,etanp,etan,                         &
                    Hprv,hv,etav,                             &
!                   --------------------------------------------
                    nRK)
          endif
!        ________________________________________
!        Hydro: Calculate velocity & pressure         
         if (RunHydro.eq.1) then		
            call hydro(                                       &
!                   -------------------------------------------
                    etanp,                                    &
                    alphafnp,ufnp,vfnp,wfnp,pfnp,             &
                    alphasnp,usnp,vsnp,wsnp,psnp,             &
!                   -------------------------------------------
                    etan,                                     &
                    alphaf,uf,vf,wf,pf,                       &
                    alphas,us,vs,ws,ps,                       &
!                   -------------------------------------------
                    etav,                                     &
                    alphafv,ufv,vfv,wfv,pfv,                  &
                    alphasv,usv,vsv,wsv,psv,                  &
!                   -------------------------------------------
                    Dpfnp,Dpf,Dpfv,                           &
!                   -------------------------------------------
                    rhof,viscof,                              &
                    rhos,viscos,                              &
!                   -------------------------------------------
                    rhofv,viscofv,                            &
                    rhosv,viscosv,                            &
!                   -------------------------------------------
                    Hpr,h,wfsurf,wssurf,                      &
!                   -------------------------------------------
                    Hprv,hv,wfsurfv,wssurfv,                  &
!                   -------------------------------------------
                    xc,yc,sig,dsig,                           &
!                   -------------------------------------------
                    xv,yv,sigv,dsigv,                         &
!                   -------------------------------------------
                    No_cp,No_vp,nbe,nbev,                     &
!                   -------------------------------------------
                    mask,                                     &
!                   -------------------------------------------
                    nRK)
         endif	
!         _______________________________________
!        |                                       |
!        |   Update RK2 values:(known)=(n+1)     |
!        |_______________________________________|

         do i=1,N_CELL
!           ----------------------------
!           Free surface  
            eta(i) = etanp(i)
            do k=1,NZ
!              -------------------------
!              Volume Fraction  
               alphaf(i,k) = alphafnp(i,k)
               alphas(i,k) = alphasnp(i,k)
!              -------------------------
!              Velocity
               uf(i,k) = ufnp(i,k)
               vf(i,k) = vfnp(i,k)
               wf(i,k) = wfnp(i,k)
               us(i,k) = usnp(i,k)
               vs(i,k) = vsnp(i,k)
               ws(i,k) = wsnp(i,k)
!              -------------------------
!              Pressure
               pf(i,k) = pfnp(i,k)
               ps(i,k) = psnp(i,k)
           enddo  
         enddo 

      ENDDO
!      ________________________________________________________ 
!     |      |                                                 |
!     | 3.2.3|  Runge-Kutta: Final solution of the RK2 formula |
!     |______|_________________________________________________|

      do i=1,N_CELL
!        ----------------------------
!        Free surface  
         etanp(i) = 0.5d0*(etan(i)+etanp(i))
         do k=1,NZ
!           -------------------------
!           Volume Fraction  
            alphafnp(i,k) = 0.5d0*(alphafn(i,k)+alphafnp(i,k))
            alphasnp(i,k) = 0.5d0*(alphasn(i,k)+alphasnp(i,k))
!           -------------------------
!           Velocity
            ufnp(i,k) = 0.5d0*(ufn(i,k)+ufnp(i,k)) 
            vfnp(i,k) = 0.5d0*(vfn(i,k)+vfnp(i,k))
            wfnp(i,k) = 0.5d0*(wfn(i,k)+wfnp(i,k))
            usnp(i,k) = 0.5d0*(usn(i,k)+usnp(i,k))
            vsnp(i,k) = 0.5d0*(vsn(i,k)+vsnp(i,k))
            wsnp(i,k) = 0.5d0*(wsn(i,k)+wsnp(i,k))
!           -------------------------
!           Pressure
            pfnp(i,k) = 0.5d0*(pfn(i,k)+pfnp(i,k))
            psnp(i,k) = 0.5d0*(psn(i,k)+psnp(i,k))
         enddo
      enddo 
!     _________________________________________________________
!     In the case of only one-step case  

      if (FinRK==1) then
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
             print*,' '
             print*,'          ---------------------------------------'
             print*,'          ------------ Using RK-1 ---------------'
             print*,'          ---------------------------------------'
             print*,' '
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
             IF (rang_topo.eq.0) THEN
             print*,' '
             print*,'          ======================================='
             print*,'          ============ Using RK-1 ==============='
             print*,'          ======================================='
             print*,' '
             ENDIF
#        endif
!        =============== END ================    
!        ====================================
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do i=1,N_CELL
            etanp(i) = eta(i)
            do k=1,NZ  
               alphafnp(i,k) = alphaf(i,k)
               ufnp(i,k) = uf(i,k) 
               vfnp(i,k) = vf(i,k)
               wfnp(i,k) = wf(i,k)
               pfnp(i,k) = pf(i,k)
!              ------------------
               alphasnp(i,k) = alphas(i,k)
               usnp(i,k) = us(i,k) 
               vsnp(i,k) = vs(i,k)
               wsnp(i,k) = ws(i,k)
               psnp(i,k) = ps(i,k)
            enddo
         enddo
      endif

      ! Remove when calculate the velocity
      do i=1,N_CELL
         etanp(i) = eta(i)
      !   do k=1,NZ  
      !      ufnp(i,k) = uf(i,k) 
      !      vfnp(i,k) = vf(i,k)
      !      wfnp(i,k) = wf(i,k)
      !  !enddo
      enddo 
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.3  |     Time loop: Update the variables & input box        |
!     |______|________________________________________________________|  

      call LoopTimeUpdate(etan,                                &
                          alphafn,ufn,vfn,wfn,pfn,             &
                          alphasn,usn,vsn,wsn,psn,             & 
                          xct,yct,zct,                         &  
                          xvt,yvt,zvt,                         &                                 
                          mask,                                & 
                          etanp,                               &
                          alphafnp,ufnp,vfnp,wfnp,pfnp,        &
                          alphasnp,usnp,vsnp,wsnp,psnp,        &
                          xc,yc,sig,Hpr,h,                     &
                          xv,yv,sigv,Hprv,hv,                  &
                          No_cp,nbe)  
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.4  |     Time loop:  Errors during simulation               |
!     |______|________________________________________________________|


      if (RunTimeExample.eq.1) then
         if (RUNTestNSEqn.eq.1) then
            !call testTimeErrorNS(ufnp,aux1,usnp,ufv,auxv1,usv,     &
            !                     vfnp,aux2,vsnp,vfv,auxv2,vsv,     &
            !                     wfnp,aux3,wsnp,wfv,auxv3,wsv,     &
            !                     pfnp,aux4,psnp,pfv,auxv4,psv,     &  
            !                     xc,yc,sig,dsig,No_cp,nbe,         &
            !                     xv,yv,sigv,dsigv,No_vp,nbev) 
            call testTimeErrorNS(ufnp,usnp,aux1,ufv,usv,auxv1,     &
                                 vfnp,vsnp,aux2,vfv,vsv,auxv2,     &
                                 wfnp,wsnp,aux3,wfv,wsv,auxv3,     &
                                 pfnp,psnp,aux4,pfv,psv,auxv4,     &  
                                 xc,yc,sig,dsig,No_cp,nbe,         &
                                 xv,yv,sigv,dsigv,No_vp,nbev) 
         else
            call testTimeError(ufnp,vfnp,wfnp,                     &
                               ufv,vfv,wfv,                        &
                               xc,yc,sig,dsig,No_cp,nbe,           &
                               xv,yv,sigv,dsigv,No_vp,nbev) 
         endif
      endif
      !________________________________________________________
      ! Free surface error
      if (RunFreeSurface.eq.1) then	
          call FS_testTimeError(ufnp,usnp,aux1,ufv,usv,auxv1,     &
                                vfnp,vsnp,aux2,vfv,vsv,auxv2,     &
                                wfnp,wsnp,aux3,wfv,wsv,auxv3,     &
                                pfnp,psnp,aux4,pfv,psv,auxv4,     &  
                                xc,yc,sig,dsig,No_cp,nbe,         &
                                xv,yv,sigv,dsigv,No_vp,nbev,      &
                                Hpr,h,etanp,                      &
                                Hprv,hv,etav,                     &
                                No_sp) 
      endif

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.5  |    Time loop:  Saving results during simulation        |
!     |______|________________________________________________________|


      if ((time.ge.tInisave).and.&
         (1d-08.ge.dtsave-(time-LastTimeSave))) then
!        ---------------------------------------------------
!        Last time saved  
         SaveCounter = SaveCounter+1
         LastTimeSave = SaveCounter*dtsave
!        ---------------------------------------------------
!        Save Tecplot file vertices 
         if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,           &
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Tecplot file cell-centers  and centers-vertex 
         if (ChooseExit.eq.2) then   
            call SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xct,yct,zct,No_cp,                 &
                               ic1tec,ic2tec,ic3tec)
!        ---------------------------------------------------
!        Save Tecplot file cell-centers & vertex 
            call SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,&
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xct,yct,zct,No_cp,                 & 
                               alphafv,ufv,vfv,wfv,pfv,           & 
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
         endif
!        -------------------------------------------------------
!        Save Re-start data  
         if (IrestartOUT.eq.1) then
            call restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp,        & 
                             alphasnp,usnp,vsnp,wsnp,psnp,        &
                             zct)
         endif
      endif

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.6  |    Time loop:  Criteria to finish time simulations     |
!     |______|________________________________________________________|   

      if ((tfin-time).le.1d-08) goto 9999
      goto 10

9999  continue


!*********************************************************************!
!                                                                     !
!                          4) Finalization                            !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.1  |                Saving final results                    |
!     |______|________________________________________________________| 

      if ((tfin-LastTimeSave).ge.1d-8) then
!        ---------------------------------------------------
!        Save Tecplot file vertices  
         if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
            print*,'   --- Tecplot final time solution printed ---'
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,           & 
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Tecplot file cell centers 
         if (ChooseExit.eq.2) then
            call SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xct,yct,zct,No_cp,                 &
                               ic1tec,ic2tec,ic3tec)
!        ---------------------------------------------------
!        Save Tecplot file cell centers & vertex 
            call SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,&
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xct,yct,zct,No_cp,                 & 
                               alphafv,ufv,vfv,wfv,pfv,           & 
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Re-start data 
         if (IrestartOUT.eq.1) then 
             call restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp,       & 
                              alphasnp,usnp,vsnp,wsnp,psnp,       &
                              zct)
         endif
      endif

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.2  |               Close time error file                    |
!     |______|________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         if (RunTimeExample.eq.1) then
            close(8100)
         endif
         if (RunFreeSurface.eq.1) then
            close(9100)
            close(9101)
            close(9102)
         endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         if (RunTimeExample.eq.1) then
            close(8100)
         endif
         if (RunFreeSurface.eq.1) then
            close(9100)
            close(9101)
            close(9102)
         endif
         ENDIF
#     endif
!     =============== END ================    
!     ====================================
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.3  |                  Final simulation time                 |
!     |______|________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         print*,'                                                        '
         print*,'    ___________________________________________________ '
         print*,'   |                                                   |'
         print*,'   |                   FINALIZATION                    |'
         print*,'   |___________________________________________________|'
         print*,'                                                        '
         print*,'     -------------------------------------------------- '
         print*,'                      Simulation Time                   '
         print*,'     -------------------------------------------------- '
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.1) THEN
         print*,'                                                        '
         print*,'   ==================================================== '
         print*,'                    MPI: FINALIZATION                   '
         print*,'   ==================================================== '
         print*,'                                                        '
         print*,'   ____________________________________________________ '
         print*,'                                                        '
         print*,'                      Simulation Time                   '
         print*,'   ____________________________________________________ '
         print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================    
!     ====================================

!     -------------------------------------------------------
!     Calculating the simulation time     
      tcpu  = etime(tt)
      idhr  = tcpu/3600
      idmin = tcpu/60-idhr*60
      idsec = tcpu-(idhr*3600+idmin*60)
!     -------------------------------------------------------
!     Print CPU Time per point    
      nt=(tfin-time0)/dt+1
      usepas=tcpu/nt/(N_CELL*NZ)*1000

!     -------------------------------------------------------
!     Display elapsed time (CPU)

916   format(t6,'PROC:',i3,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
917   format(t6,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
918   format(t6,'          user : ',f10.3,' sec')
919   format(t6,'        system : ',f10.3,' sec')  
920   format(t6,'Time per point : ',f10.3,' msec CPU /pas/point ')

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(*,917) tcpu,idhr,idmin,idsec
         write(*,918) tt(1)
         write(*,919) tt(2)
         print*,'  '
         write(*,920) usepas
         print*,'                                                        '
         print*,'     -------------------------------------------------- '
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         write(*,916) rang_topo,tcpu,idhr,idmin,idsec
         call MPI_ALLREDUCE(tcpu,MAXtcpu,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                          comm3D,code)
         IF (rang_topo.eq.0) THEN
            print*,' '
            print*,'    Maximum time = ',MAXtcpu 
            print*,' '
            print*,'   ==================================================== '
            print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================    
!     ====================================
   
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.4  |             Free memory: deallocating                  |
!     |______|________________________________________________________|

      call dealloc_variables
      call dealloc_geometry 
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.5  |                 Closing programns                      |
!     |______|________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call finalisation_mpi
#     endif
!     =============== END ================    
!     ====================================    

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'------------------------------------------------------------' 
      print*,'                      ***  END  ***                         '
      print*,'------------------------------------------------------------'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
      IF (rang_topo.eq.0) THEN
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'============================================================' 
      print*,'                  ***  END MPI PROGRAM ***                  '
      print*,'============================================================'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
      ENDIF
#     endif
!     =============== END ================ 

      END PROGRAM 

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                            END OF nsmp3D                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
