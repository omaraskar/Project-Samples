!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          TESTING SUBROUTINE                         !
!              Advection, Diffusion and Poisson components            !
!                             May 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

       SUBROUTINE testTimeExamples(ufnp,ufn,uf,ufv,                   &
                                   vfnp,vfn,vf,vfv,                   &
                                   wfnp,wfn,wf,wfv,                   &
                                   pfnp,pfn,pf,pfv,                   &
!                                  ------------------------------------
                                   usnp,usv,                          &
                                   vsnp,vsv,                          &
                                   wsnp,wsv,                          &
                                   psnp,psv,                          &
!                                  ------------------------------------
                                   xc,yc,sig,dsig,No_cp,nbe,          &
!                                  ------------------------------------
                                   xv,yv,sigv,dsigv,No_vp,nbev,       &
!                                  ------------------------------------      
                                   Hpr,h,etan,                        &
!                                  ------------------------------------
                                   Hprv,hv,etav,                      &
!                                  ------------------------------------
                                   nRK)

!---------------------------------------------------------------------!   
!                                                                     !
!     This subroutine tests several subroutines developed to solve    !
!     the Navier-Stokes equations.                                    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !  
!  |____________|____________|_____________________________________|  ! 
!  | <-- ufnp   |(N_CELL,NZ) | Solution of the test at t(n+1)      |  !
!  | <-- ufn    |(N_CELL,NZ) | Solution of the test at t(n)        |  !      
!  | <-- uf     |(N_CELL,NZ) | Solution of the test at RK step     |  !  
!  | <-- ufv    |(N_VERT,NZ) | Solution of the test at the vertices|  ! 
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbev   |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!  | --> nRK    | integer    | Rouge-Kutta iteration               |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - testTimeAdvEqn      (TestExamples/testTimeAdvEqn.F90    ) |  !
!  |   - testTimeDiffEqn     (TestExamples/testTimeDiffEqn.F90   ) |  !
!  |   - testTimeAdvDiffEqn  (TestExamples/testTimeAdvDiffEqn.F90) |  !
!  |   - testTimeNSEqn       (TestExamples/testTimeNSEqn.F90     ) |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!   ----  Parameters                                                  !
!        Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            implicit none
#     endif
!     =============== END ================    
!     ====================================  
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|
                                       
      real*8, dimension(:,:):: ufnp(N_CELL,NZ) 
      real*8, dimension(:,:):: ufn(N_CELL,NZ)                 
      real*8, dimension(:,:):: uf(N_CELL,NZ)
      real*8, dimension(:,:):: ufv(N_VERT,NZ-1)
      real*8, dimension(:,:):: vfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: vfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: vf(N_CELL,NZ)
      real*8, dimension(:,:):: vfv(N_VERT,NZ-1)
      real*8, dimension(:,:):: wfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: wfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: wf(N_CELL,NZ)
      real*8, dimension(:,:):: wfv(N_VERT,NZ-1)
      real*8, dimension(:,:):: pfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: pfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: pf(N_CELL,NZ)
      real*8, dimension(:,:):: pfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:,:):: usnp(N_CELL,NZ) 
      real*8, dimension(:,:):: usv(N_VERT,NZ-1)
      real*8, dimension(:,:):: vsnp(N_CELL,NZ) 
      real*8, dimension(:,:):: vsv(N_VERT,NZ-1)
      real*8, dimension(:,:):: wsnp(N_CELL,NZ) 
      real*8, dimension(:,:):: wsv(N_VERT,NZ-1)
      real*8, dimension(:,:):: psnp(N_CELL,NZ) 
      real*8, dimension(:,:):: psv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     --------------------------------------
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      integer :: nRK 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
              
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gam(N_CELL,NZ)
!     --------------------------------------
      real*8,dimension(:) :: dHprdx(N_CELL)
      real*8,dimension(:) :: dHprdy(N_CELL)
      real*8,dimension(:) :: dhdx(N_CELL)
      real*8,dimension(:) :: dhdy(N_CELL)
!     --------------------------------------
      real*8 :: funEta,funh
!     --------------------------------------
      real*8 :: x,y,z
      integer:: tagU

!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
          write(*,'(t8,60a)'),'----> Begin subroutine: testing time examples'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Sigma transform for the examples                   !
!                                                                     !
!*********************************************************************!

!     ________________________________________
!     Eta  
      do i=1,N_CELL 
         etan(i)  = funEta(xc(i),yc(i))
      enddo
      do nv=1,N_VERT 
         etav(nv) = funEta(xv(nv),yv(nv)) 
      enddo
!     ________________________________________
!     Update depth: H = h + eta 
!     -------------------------
!     cell center 
      do i=1,N_CELL 
         Hpr(i) = h(i) + etan(i) 
      enddo
!     -------------------------
!     vertex
      do nv=1,N_VERT 
         Hprv(nv) = hv(nv) + etav(nv) 
      enddo
!     ________________________________________
!     Gradients
      call grandientLSM2D(dHprdx,dHprdy, &
                          Hpr,No_cp,nbe) 
      call grandientLSM2D(dhdx,dhdy, &
                          h,No_cp,nbe) 
!     ________________________________________
!     Sigma values
      do k=1,NZ
         do i=1,N_CELL 
            sigmax(i,k) = 1.0d0/Hpr(i)*(dhdx(i)-sig(k)*dHprdx(i)) 
            sigmay(i,k) = 1.0d0/Hpr(i)*(dhdy(i)-sig(k)*dHprdy(i)) 
            sigmaz(i,k) = 1.0d0/Hpr(i)                            
         enddo
      enddo

!     ________________________________________
!     Box domain case
      if (ChooseDomBox.eq.1) then
         do i=1,N_CELL 
            h(i) = 0.0d0
            Hpr(i) = 1.0d0 
         enddo
         do nv=1,N_VERT 
            hv(nv)   = 0.0d0
            Hprv(nv) = 1.0d0 
         enddo
         do k=1,NZ
            do i=1,N_CELL 
               sigmax(i,k) = 0.0d0
               sigmay(i,k) = 0.0d0
               sigmaz(i,k) = 1.0d0  
            enddo
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                       >>>>> TEST ZONE <<<<<                         !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                 Testing Time Equations                 |
!     |________________________________________________________|

!     __________________________________________________________
!     >>>> TEST: RK-2 time approximation

      if (RUNtestRK2approx.eq.1)  then
         do i=1,N_CELL0
            do k=1,NZ
               ufnp(i,k) = uf(i,k) + dt*dcos(time+dt*(nRK-2))
            enddo
         enddo 
      endif   
!     __________________________________________________________
!     >>>> TEST: Advection equation

      if (RUNtestAdvEqn.eq.1)  then
         call testTimeAdvEqn(ufnp,ufv,                        &  
                             uf,xc,yc,sig,dsig,No_cp,nbe,     &
                             ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                             nRK)    
      endif
!     __________________________________________________________
!     >>>> TEST: Diffusion equation

      if (RUNtestDiffEqn.eq.1)  then
         call testTimeDiffEqn(ufnp,ufv,                        &  
                              uf,xc,yc,sig,dsig,No_cp,nbe,     &
                              ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                              nRK)    
      endif
!     __________________________________________________________
!     >>>> TEST: Advection-Diffusion equation

      if (RUNtestAdvDiffEqn.eq.1)  then
         call testTimeAdvDiffEqn(ufnp,ufv,                        &  
                                 uf,xc,yc,sig,dsig,No_cp,nbe,     &
                                 ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 nRK)    
      endif

!     __________________________________________________________
!     >>>> TEST: Navier-Stokes equation 3D

      if (RUNTestNSEqn.eq.1)  then
         call testTimeNSEqn(ufnp,vfnp,wfnp,pfnp,             &
                            ufv,vfv,wfv,pfv,                 & 
                            uf,vf,wf,pf,                     &
                            xc,yc,sig,dsig,No_cp,nbe,        &
                            xv,yv,sigv,dsigv,No_vp,nbev,     &
                            Hpr,h,etan,                      &
                            Hprv,hv,etav,                    &
                            nRK)    
      endif

!*********************************************************************!
!                                                                     !
!                    Finalization of the subroutine                   !
!                                                                     !
!*********************************************************************!  

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
          write(*,'(t8,60a)'),'---->   End subroutine: testing'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	 END testing Examples                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!  
