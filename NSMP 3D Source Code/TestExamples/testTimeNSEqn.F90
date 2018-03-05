!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    TEST NAVIER-STOKES EQUATION                      !
!                              Nov 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeNSEqn(ufnp,vfnp,wfnp,pfnp,                  &  
                               ufv,vfv,wfv,pfv,                      &
                               uf,vf,wf,pf,                          &
                               xc,yc,sig,dsig,No_cp,nbe,             &
                               xv,yv,sigv,dsigv,No_vp,nbev,          &
                               Hpr,h,etan,                           &
                               Hprv,hv,etav,                         &
                               nRK)               

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the 3D Navier-Stokes equations with velo-   !
!    city variables u, v, w and pressure p. We have two methods to    !
!    calculate teh pressure:                                          !
!                                                                     !
!    PressureMethod = 1 :  Poisson eqn using the projection method.   !
!    PressureMethod = 2 :  Poisson eqn using a direct method.         !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- ufnp   |(N_CELL,NZ)  | u Approximate solution cell-center|  !   
!  | <--> ufv    |(N_VERT,NZ-1)| u Approximate solution vertex     |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- vfnp   |(N_CELL,NZ)  | v Approximate solution cell-center|  !  
!  | <--> vfv    |(N_VERT,NZ-1)| v Approximate solution vertex     |  ! 
!  |_____________|_____________|___________________________________|  !
!  | <--- wfnp   |(N_CELL,NZ)  | w Approximate solution cell-center|  !   
!  | <--> wfv    |(N_VERT,NZ-1)| w Approximate solution vertex     |  ! 
!  |_____________|_____________|___________________________________|  !
!  | <--- pfnp   |(N_CELL,NZ)  | p Approximate solution cell-center|  !  
!  | <--> pfv    |(N_VERT,NZ-1)| p Approximate solution vertex     |  !  
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> uf    |(N_CELL,NZ)  | u Old solution of the equation     |  !
!  | ---> vf    |(N_CELL,NZ)  | v Old solution of the equation     |  !  
!  | ---> wf    |(N_CELL,NZ)  | w Old solution of the equation     |  !  
!  | ---> pf    |(N_CELL,NZ)  | p Old solution of the equation     |  !
!  |____________|_____________|____________________________________|  !  
!  | ---> xc,yc |(N_CELL)     | Coordinates of the cell centers    |  !
!  | ---> sig   |(NZ)         | Sigma value at the cell centers    |  !
!  | ---> dsig  |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | ---> No_cp |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | ---> nbe   |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | ---> sigv  |(NZ-1)       | sigma of the vertex points         |  !
!  | ---> dsigv |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !  
!  | ---> No_vp |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | ---> nbev  |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - AdvDiffEqn3D                ( AdvDiffEqn3D.F90 )          |  !
!  |   - Poisson                     ( Poisson.F90 )               |  !
!  |   - grandientLSM                ( grandientLSM.F90 )          |  !
!  |   - interpolation3D             ( interpolation3D.F90 )       |  !
!  |   - BCvelcenter3D               ( BCvelocity.F90 )            |  !
!  |   - BCvelvertex3D               ( BCvelocity.F90 )            |  !
!  |   - BCprescenter3D              ( BCpressure.F90 )            |  !
!  |   - BCpresvertex3D              ( BCpressure.F90 )            |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!     ____________________________________
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|

      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: uf(N_CELL,NZ)
      real*8,dimension(:,:) :: vf(N_CELL,NZ)
      real*8,dimension(:,:) :: wf(N_CELL,NZ)
      real*8,dimension(:,:) :: pf(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

!     --------------------------------------
      real*8,dimension(:,:),allocatable :: ufstar,vfstar,wfstar
      real*8,dimension(:,:),allocatable :: Dpfnp,Dpf,Dpfnpv,Dpfv
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: dudx,dudy,dudz
      real*8,dimension(:,:),allocatable :: dvdx,dvdy,dvdz
      real*8,dimension(:,:),allocatable :: dwdx,dwdy,dwdz
      real*8,dimension(:,:),allocatable :: dpdx,dpdy,dpdz
      real*8,dimension(:,:),allocatable :: dDpdx,dDpdy,dDpdz
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: Gamu,Gamv,Gamw,Gamp
      real*8,dimension(:,:),allocatable :: rhsu,rhsv,rhsw,rhsp
!     -------------------------------------
      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phiE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8,dimension(:) :: SaveErrorMax(4)
      real*8,dimension(:) :: SaveErrorSum(4)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
!     --------------------------------------     
      real*8 :: x,y,z,uu,vv,ww
      real*8 :: funExamNSu,funExamNSv,funExamNSw
      real*8 :: funExamNSp,funExamNSrhsp 
      integer:: tagBC,UseThis,DisplayThis,s
!     --------------------------------------
      integer, parameter :: PressureMethod  = 2
      integer, parameter :: TestOnlyAdvDiff = 0
!     --------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TEST NS Equation'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(ufstar(N_CELL,NZ),vfstar(N_CELL,NZ),wfstar(N_CELL,NZ), &
               Dpfnp(N_CELL,NZ),Dpf(N_CELL,NZ),                       &
               Dpfnpv(N_VERT,NZ-1),Dpfv(N_VERT,NZ-1),                 &            
               dudx(N_CELL,NZ),dudy(N_CELL,NZ),dudz(N_CELL,NZ),       &
               dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvdz(N_CELL,NZ),       &
               dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwdz(N_CELL,NZ),       &
               dpdx(N_CELL,NZ),dpdy(N_CELL,NZ),dpdz(N_CELL,NZ),       &
               dDpdx(N_CELL,NZ),dDpdy(N_CELL,NZ),dDpdz(N_CELL,NZ),    &
               Gamu(N_CELL,NZ),Gamv(N_CELL,NZ),Gamw(N_CELL,NZ),       &
               Gamp(N_CELL,NZ),                                       &
               rhsu(N_CELL,NZ),rhsv(N_CELL,NZ),rhsw(N_CELL,NZ),       &
               rhsp(N_CELL,NZ))

!      ________________________________________________________
!     |                                                        |
!     |             Total Depth & sigma transformation         |
!     |________________________________________________________|

      if (ChooseDomBox.eq.1) then
         do i=1,N_CELL 
            Hpr(i) = 1.0d0 
            do k=1,NZ
               sigmax(i,k) = 0.0d0
               sigmay(i,k) = 0.0d0
               sigmaz(i,k) = 1.0d0  
            enddo
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                       Navier-Stokes solution                        !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |          ----  3) PREDICTION STEP: fluid  ----         |
!     |                                                        |        
!     |                    °  uf^(*)                           | 
!     |                    °  vf^(*)                           |   
!     |                    °  wf^(*)                           |
!     |________________________________________________________|
!      ********************************************************
!      _______________________________________________
!     |                                               |
!     |  3.0)  Coefficients: diffusive coeff = Gam    |
!     |                      righ-hand side  = rhs    |
!     |_______________________________________________|

!     _________________________________________________
!     Diffusive Coefficients: Gamu, Gamv, Gamw
      do k=1,NZ
         do i=1,N_CELL  
            Gamu(i,k) = 1.0d0/Re  
            Gamv(i,k) = 1.0d0/Re  
            Gamw(i,k) = 1.0d0/Re  
         enddo
      enddo
!     _________________________________________________
!     Right-hand side: rhs
!     ----------------------------- 
      call grandientLSM(dpdx,dpdy,dpdz,pf,No_cp,nbe,sig) 
!     ----------------------------- 
      do k=1,NZ
         do i=1,N_CELL  
            rhsu(i,k) = -dpdx(i,k)
            rhsv(i,k) = -dpdy(i,k)
            rhsw(i,k) = -dpdz(i,k) 
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  3.1) Update velocity components: uf,vf,wf^(*)|
!     |_______________________________________________|

!     ----------------------------- 
      call AdvDiffVelocity(ufstar,ufv,                      & 
                           rhsu,Gamu,Gamu,Gamu,uf,vf,wf,pf, &  
                           uf,xc,yc,sig,dsig,No_cp,nbe,     &
                           ufv,xv,yv,sigv,dsigv,No_vp,nbev,1)
!     ----------------------------- 
      call AdvDiffVelocity(vfstar,vfv,                      & 
                           rhsv,Gamv,Gamv,Gamv,uf,vf,wf,pf, & 
                           vf,xc,yc,sig,dsig,No_cp,nbe,     &
                           vfv,xv,yv,sigv,dsigv,No_vp,nbev,2)
!     ----------------------------- 
      call AdvDiffVelocity(wfstar,wfv,                      & 
                           rhsw,Gamw,Gamw,Gamw,uf,vf,wf,pf, & 
                           wf,xc,yc,sig,dsig,No_cp,nbe,     &
                           wfv,xv,yv,sigv,dsigv,No_vp,nbev,3)

!     _________________________________________________
!     -----    OPTIONAL TEST: Exact velocity     ------
      UseThis = 0
      if (UseThis.eq.1) then
         do k=1,NZ
            do i=1,N_CELL
               x = xc(i)
               y = yc(i)
               z = sig(k)
               ufstar(i,k) = funExamNSu(x,y,z,time)
               vfstar(i,k) = funExamNSv(x,y,z,time)
               wfstar(i,k) = funExamNSw(x,y,z,time)
            enddo
         enddo
      endif
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |          ----  4) PROJECTION STEP: fluid   ----        |
!     |                                                        |        
!     |                    °  pf^(new)                         |
!     |________________________________________________________|
!      ********************************************************
!      _______________________________________________
!     |                                               |
!     |  4.0)  Coefficients: diffusive coeff = Gamp   |
!     |                      righ-hand side  = rhsp   |
!     |_______________________________________________|

!     __________________________________________________
!     Diffusive Coefficients: Gamp 
      do k=1,NZ
         do i=1,N_CELL  
            Gamp(i,k) = 1.0d0
         enddo
      enddo
!     _________________________________________________
!     Right-hand side 
!        ----------------------------- 
      call grandientLSM(dudx,dudy,dudz,ufstar,No_cp,nbe,sig) 
      call grandientLSM(dvdx,dvdy,dvdz,vfstar,No_cp,nbe,sig) 
      call grandientLSM(dwdx,dwdy,dwdz,wfstar,No_cp,nbe,sig) 
!        ----------------------------- 
      IF (PressureMethod.eq.1) THEN
         do k=1,NZ
            do i=1,N_CELL  
               rhsp(i,k) = (dudx(i,k)+dvdy(i,k)+dwdz(i,k))/dt
            enddo
         enddo
      ELSEIF (PressureMethod.eq.2) THEN
         do k=1,NZ
            do i=1,N_CELL  
               rhsp(i,k) = -2.0d0*                                 &
                           (dudy(i,k)*dvdx(i,k)-dudx(i,k)*dvdy(i,k)&
                           +dudz(i,k)*dwdx(i,k)-dudx(i,k)*dwdz(i,k)&
                           +dvdz(i,k)*dwdy(i,k)-dvdy(i,k)*dwdz(i,k)) 
            enddo
         enddo
      ENDIF
!     _________________________________________________
!     -----      OPTIONAL TEST: Exact rhsp       ------
      UseThis = 0
      if (UseThis.eq.1) then
         do k=1,NZ
            do i=1,N_CELL  
               x  = xc(i)
               y  = yc(i)
               z  = sig(k)
               uu = ufstar(i,k)
               vv = vfstar(i,k)
               ww = wfstar(i,k) 
               rhsp(i,k) = funExamNSrhsp(x,y,z,uu,vv,ww,time)
            enddo
         enddo
      endif
!      _______________________________________________
!     |                                               |
!     |  4.1)    Poisson solver: Div(Dp)=rhs          |
!     |            where Dp = pf^(n+1)-pf             |
!     |_______________________________________________|

      IF (PressureMethod.eq.1) THEN
         call Poisson(Dpfnp,Dpfv,                    &
                      rhsp,Gamp,Gamp,Gamp,           & 
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,   &
                      Hpr,h,etan,                    &
                      Hprv,hv,etav) 

      ELSEIF (PressureMethod.eq.2) THEN
         call Poisson(pfnp,pfv,                      &
                      rhsp,Gamp,Gamp,Gamp,           & 
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,   &
                      Hpr,h,etan,                    &
                      Hprv,hv,etav) 
      ENDIF

!      _______________________________________________
!     |                                               |
!     |  4.2)  Pressure solution: pf^(n+1) = Dp + pf  |
!     |_______________________________________________|

      IF (PressureMethod.eq.1) THEN
!        ______________________________________________
!        Update cell-centers
         do k=1,NZ
            do i=1,N_CELL  
               pfnp(i,k) = Dpfnp(i,k) + pf(i,k) 
            enddo
         enddo
!        ______________________________________________
!        Update vertex values
         do k=1,NZ-1
            do nv=1,N_VERT  
               pfv(nv,k) = Dpfv(nv,k) + pfv(nv,k) 
            enddo
         enddo
       ENDIF

!     _________________________________________________
!     -----    OPTIONAL TEST: Exact pressure p   ------
      UseThis = 0
      if (UseThis.eq.1) then
         do k=1,NZ
            do i=1,N_CELL  
               x  = xc(i)
               y  = yc(i)
               z  = sig(k)
               pfnp(i,k) = funExamNSp(x,y,z,time)
            enddo
         enddo
      endif
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |           ---- 5) CORRECTION STEP: fluid ----          |
!     |                                                        | 
!     |                    °  uf^(new)                         | 
!     |                    °  vf^(new)                         |   
!     |                    °  wf^(new)                         |
!     |________________________________________________________|
!      ********************************************************
!      _______________________________________________
!     |                                               |
!     |  5.2)  Update velocity components             |
!     |_______________________________________________|

      IF (PressureMethod.eq.1) THEN
!        ______________________________________________
!        Update velocity: cell-centers
!        ----------------------------- 
         call grandientLSM(dDpdx,dDpdy,dDpdz,Dpfnp,No_cp,nbe,sig) 
!        ----------------------------- 
         do k=1,NZ
            do i=1,N_CELL0
               ufnp(i,k) = ufstar(i,k)-dt*dDpdx(i,k)
               vfnp(i,k) = vfstar(i,k)-dt*dDpdy(i,k)
               wfnp(i,k) = wfstar(i,k)-dt*dDpdz(i,k) 
            enddo
         enddo
      ELSEIF (PressureMethod.eq.2) THEN
!        ______________________________________________
!        Right-hand side: rhs
!        ----------------------------- 
         call grandientLSM(dpdx,dpdy,dpdz,pf,No_cp,nbe,sig) 
!        ----------------------------- 
         do k=1,NZ
            do i=1,N_CELL  
               rhsu(i,k) = -dpdx(i,k)
               rhsv(i,k) = -dpdy(i,k)
               rhsw(i,k) = -dpdz(i,k) 
            enddo
         enddo
!        ______________________________________________
!        Update velocity components
!        ----------------------------- 
         call AdvDiffVelocity(ufnp,ufv,                         & 
                              rhsu,Gamu,Gamu,Gamu,uf,vf,wf,pfnp,&  
                              uf,xc,yc,sig,dsig,No_cp,nbe,      &
                              ufv,xv,yv,sigv,dsigv,No_vp,nbev,1)
!        -----------------------------
         call AdvDiffVelocity(vfnp,vfv,                         & 
                              rhsv,Gamv,Gamv,Gamv,uf,vf,wf,pfnp,&  
                              vf,xc,yc,sig,dsig,No_cp,nbe,      &
                              vfv,xv,yv,sigv,dsigv,No_vp,nbev,2)
!        -----------------------------
         call AdvDiffvelocity(wfnp,wfv,                         & 
                              rhsw,Gamw,Gamw,Gamw,uf,vf,wf,pfnp,&
                              wf,xc,yc,sig,dsig,No_cp,nbe,      &
                              wfv,xv,yv,sigv,dsigv,No_vp,nbev,3)
      ENDIF

!*********************************************************************!
!                                                                     !
!                     Test Adv-Diff velocity (Exact p)                !
!                                                                     !
!*********************************************************************!

     IF (TestOnlyAdvDiff.eq.1) THEN
!        _________________________________________________
!        Diffusive Coefficients: Gamu, Gamv, Gamw
         do k=1,NZ
            do i=1,N_CELL  
               Gamu(i,k) = 1.0d0/Re  
               Gamv(i,k) = 1.0d0/Re  
               Gamw(i,k) = 1.0d0/Re  
            enddo
         enddo
!        _________________________________________________
!        Exact solutions
         do k=1,NZ
            do i=1,N_CELL  
               x  = xc(i)
               y  = yc(i)
               z  = sig(k)
               pfnp(i,k)= funExamNSp(x,y,z,time)
            enddo
         enddo
!        ______________________________________________
!        Gradient of p
         call grandientLSM(dpdx,dpdy,dpdz,pf,No_cp,nbe,sig) 
!        ______________________________________________
!        Right-hand side: rhs
         do k=1,NZ
            do i=1,N_CELL  
               rhsu(i,k) = -dpdx(i,k)
               rhsv(i,k) = -dpdy(i,k)
               rhsw(i,k) = -dpdz(i,k) 
            enddo
         enddo
!        ______________________________________________
!        Update velocity components
         call AdvDiffVelocity(ufnp,ufv,                         & 
                              rhsu,Gamu,Gamu,Gamu,uf,vf,wf,pfnp,&  
                              uf,xc,yc,sig,dsig,No_cp,nbe,      &
                              ufv,xv,yv,sigv,dsigv,No_vp,nbev,1)
!        -----------------------------
         call AdvDiffVelocity(vfnp,vfv,                         & 
                              rhsv,Gamv,Gamv,Gamv,uf,vf,wf,pfnp,&  
                              vf,xc,yc,sig,dsig,No_cp,nbe,      &
                              vfv,xv,yv,sigv,dsigv,No_vp,nbev,2)
!        -----------------------------
         call AdvDiffvelocity(wfnp,wfv,                         & 
                              rhsw,Gamw,Gamw,Gamw,uf,vf,wf,pfnp,&
                              wf,xc,yc,sig,dsig,No_cp,nbe,      &
                              wfv,xv,yv,sigv,dsigv,No_vp,nbev,3)
      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(ufstar,vfstar,wfstar,  &
                 Dpfnp,Dpf,Dpfnpv,Dpfv, &            
                 dudx,dudy,dudz,        &
                 dvdx,dvdy,dvdz,        &
                 dwdx,dwdy,dwdz,        &
                 dpdx,dpdy,dpdz,        &
                 dDpdx,dDpdy,dDpdz,     &
                 Gamu,Gamv,Gamw,Gamp,   &
                 rhsu,rhsv,rhsw,rhsp)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TEST NS equation'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   END OF TEST NAVIER-STOKES EQUATION                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
