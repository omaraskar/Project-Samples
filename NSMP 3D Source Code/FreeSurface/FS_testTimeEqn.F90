!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                FREE SURFACE TEST NAVIER-STOKES EQUATION             !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_testTimeEqn(ufnp,ufn,uf,ufv,                   &
                                vfnp,vfn,vf,vfv,                   &
                                wfnp,wfn,wf,wfv,                   &
                                pfnp,pfn,pf,pfv,                   &
!                               ------------------------------------
                                rhof,viscof,                       &
                                rhofv,viscofv,                     &
!                               ------------------------------------
                                xc,yc,sig,dsig,No_cp,nbe,          &
                                xv,yv,sigv,dsigv,No_vp,nbev,       &
!                               ------------------------------------
                                Hpr,h,eta,etan,                    &
                                Hprv,hv,etav,                      &
!                               ------------------------------------
                                nRK)               

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the 3D Navier-Stokes equations with velo-   !
!    city variables (u,v,w) and pressure p. We have two methods to    !
!    calculate the pressure:                                          !
!                                                                     !
!    PressureMethod = 1 : Poisson eqn using the projection method.    !
!    PressureMethod = 2 : Poisson eqn using a direct method (better). !
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
!  | ---> uf    |(N_CELL,NZ)  | u solution at the current RK step  |  !
!  | ---> vf    |(N_CELL,NZ)  | v solution at the current RK step  |  !
!  | ---> wf    |(N_CELL,NZ)  | w solution at the current RK step  |  !
!  | ---> pf    |(N_CELL,NZ)  | p solution at the current RK step  |  !
!  |____________|_____________|____________________________________|  !
!  | ---> ufn   |(N_CELL,NZ)  | u solution at time n               |  !
!  | ---> vfn   |(N_CELL,NZ)  | v solution at time n               |  !
!  | ---> wfn   |(N_CELL,NZ)  | w solution at time n               |  ! 
!  | ---> pfn   |(N_CELL,NZ)  | p solution at time n               |  !
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

      real*8, dimension(:,:):: ufnp(N_CELL,NZ) 
      real*8, dimension(:,:):: ufn(N_CELL,NZ)                 
      real*8, dimension(:,:):: uf(N_CELL,NZ)
      real*8, dimension(:,:):: ufv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: vfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: vfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: vf(N_CELL,NZ)
      real*8, dimension(:,:):: vfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: wfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: wfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: wf(N_CELL,NZ)
      real*8, dimension(:,:):: wfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: pfnp(N_CELL,NZ) 
      real*8, dimension(:,:):: pfn(N_CELL,NZ)                 
      real*8, dimension(:,:):: pf(N_CELL,NZ)
      real*8, dimension(:,:):: pfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: rhof(N_CELL,NZ)
      real*8,dimension(:,:) :: viscof(N_CELL,NZ) 
!     -------------------------------------
      real*8,dimension(:,:) :: rhofv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: viscofv(N_VERT,NZ-1) 
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
      real*8, dimension(:)  :: eta(N_CELL)
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
      real*8,dimension(:,:),allocatable :: dudx,dudy,duds
      real*8,dimension(:,:),allocatable :: dvdx,dvdy,dvds
      real*8,dimension(:,:),allocatable :: dwdx,dwdy,dwds
      real*8,dimension(:,:),allocatable :: dpdx,dpdy,dpds
      real*8,dimension(:,:),allocatable :: dDpdx,dDpdy,dDpds
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: Gamp
      real*8,dimension(:,:),allocatable :: Gamux,Gamuy,Gamuz
      real*8,dimension(:,:),allocatable :: Gamvx,Gamvy,Gamvz
      real*8,dimension(:,:),allocatable :: Gamwx,Gamwy,Gamwz
      real*8,dimension(:,:),allocatable :: rhsu,rhsv,rhsw,rhsp
      real*8,dimension(:,:),allocatable :: ph,omegaf
      real*8,dimension(:,:),allocatable :: Huf,Hvf,Hwf
      real*8,dimension(:,:),allocatable :: Hufv,Hvfv,Hwfv
!     --------------------------------------
      real*8, dimension(:)  :: Hold(N_CELL)
      real*8, dimension(:)  :: Hvold(N_VERT)
!     --------------------------------------
      real*8,dimension(:) :: dHprdx(N_CELL)
      real*8,dimension(:) :: dHprdy(N_CELL)
      real*8,dimension(:) :: dhdx(N_CELL)
      real*8,dimension(:) :: dhdy(N_CELL)
      real*8,dimension(:) :: detadx(N_CELL)
      real*8,dimension(:) :: detady(N_CELL)
      real*8,dimension(:) :: detadt(N_CELL)
      real*8,dimension(:) :: d2etadx2(N_CELL)
      real*8,dimension(:) :: d2etady2(N_CELL)
      real*8,dimension(:) :: d2etadxdy(N_CELL)
      real*8,dimension(:) :: d2etadydx(N_CELL)
      real*8,dimension(:) :: dwfdtB(N_CELL)
      real*8,dimension(:) :: dwfdtT(N_CELL)
      real*8,dimension(:) :: dwfvdtB(N_VERT)
      real*8,dimension(:) :: dwfvdtT(N_VERT)
!     --------------------------------------
      real*8,dimension(:) :: etaExact(N_CELL)
      real*8,dimension(:) :: HprExact(N_CELL)
!     --------------------------------------
      real*8 :: Sdudx,Sdudy,Sdudz
      real*8 :: Sdvdx,Sdvdy,Sdvdz
      real*8 :: Sdwdx,Sdwdy,Sdwdz
      real*8 :: Sdpdx,Sdpdy,Sdpdz
      real*8 :: SdDpdx,SdDpdy,SdDpdz
!     --------------------------------------
      real*8 :: x,y,z,uu,vv,ww
      real*8 :: FS_funu,FS_funv,FS_funw,FS_funeta,etaE
      real*8 :: Cs,DD,nut,Sij,c1,c2
      real*8 :: uT,uB,vT,vB
      real*8 :: rhsph
      integer:: UseThis,DisplayThis,s,jc
      integer:: tagBCu,tagBCv,tagBCw
!     --------------------------------------
      integer, parameter :: DynamicPressure = 1
      integer, parameter :: PressureMethod  = 2
!     --------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TEST Free Surface'
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
               dudx(N_CELL,NZ),dudy(N_CELL,NZ),duds(N_CELL,NZ),       &
               dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvds(N_CELL,NZ),       &
               dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwds(N_CELL,NZ),       &
               dpdx(N_CELL,NZ),dpdy(N_CELL,NZ),dpds(N_CELL,NZ),       &
               dDpdx(N_CELL,NZ),dDpdy(N_CELL,NZ),dDpds(N_CELL,NZ),    &
               Gamux(N_CELL,NZ),Gamuy(N_CELL,NZ),Gamuz(N_CELL,NZ),    &
               Gamvx(N_CELL,NZ),Gamvy(N_CELL,NZ),Gamvz(N_CELL,NZ),    &
               Gamwx(N_CELL,NZ),Gamwy(N_CELL,NZ),Gamwz(N_CELL,NZ),    &
               Gamp(N_CELL,NZ),                                       &
               rhsu(N_CELL,NZ),rhsv(N_CELL,NZ),rhsw(N_CELL,NZ),       &
               rhsp(N_CELL,NZ),ph(N_CELL,NZ),omegaf(N_CELL,NZ),       &
               Huf(N_CELL,NZ),Hvf(N_CELL,NZ),Hwf(N_CELL,NZ),          &
               Hufv(N_VERT,NZ-1),Hvfv(N_VERT,NZ-1),Hwfv(N_VERT,NZ-1))

!*********************************************************************!
!                                                                     !
!                        Free-surface solution                        !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |           ---- 1) UPDATE FREE SURFACE    ----          |
!     |                     °  (eta)^(new)                     |        
!     |                     °  (Hpr)^(new)                     |
!     |________________________________________________________|
!      ********************************************************

!      _______________________________________________
!     |                                               |
!     |  1.1)  Update H using free surface equation   |
!     |_______________________________________________|

      call FreeSurface(Hpr,h,eta,                   &
                       Hprv,hv,etav,                &
                       uf,vf,                       & 
                       xc,yc,sig,dsig,No_cp,nbe,    &
                       xv,yv,sigv,dsigv,No_vp,nbev, &
                       nRK)   

!     ________________________________________
!     Exact eta and Hpr 
      do i=1,N_CELL 
         x = xc(i)
         y = yc(i)
         etaExact(i) = FS_funeta(x,y,time)
         HprExact(i) = etaExact(i) + h(i)
      enddo

!      _______________________________________________
!     |                                               |
!     |  1.2)  Sigma transformation                   |
!     |_______________________________________________|

!     ________________________________________
!     Derivatives 
      do i=1,N_CELL 
         detadt(i) = (eta(i)-etan(i))/dt
      enddo
      call grandientLSM2D(dHprdx,dHprdy,Hpr,No_cp,nbe) 
      call grandientLSM2D(detadx,detady,eta,No_cp,nbe) 
      call grandientLSM2D(dhdx,dhdy,h,No_cp,nbe) 
      call grandientLSM2D(d2etadx2,d2etadxdy,detadx,No_cp,nbe) 
      call grandientLSM2D(d2etadydx,d2etady2,detady,No_cp,nbe) 
!     ________________________________________
!     Sigma transform
      do k=1,NZ
         do i=1,N_CELL 
            sigmat(i,k) = 1.0d0/Hpr(i)*(-sig(k)*detadt(i))
            sigmax(i,k) = 1.0d0/Hpr(i)*(dhdx(i)-sig(k)*dHprdx(i)) 
            sigmay(i,k) = 1.0d0/Hpr(i)*(dhdy(i)-sig(k)*dHprdy(i)) 
            sigmaz(i,k) = 1.0d0/Hpr(i)                           
         enddo
      enddo
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |       ----  2) VELOCITY: PREDICTION STEP  ----         |       
!     |                    °  uf^(*)                           | 
!     |                    °  vf^(*)                           |   
!     |                    °  wf^(*)                           |
!     |________________________________________________________|
!      ********************************************************
!      _______________________________________________
!     |                                               |
!     |  2.1)  Coefficients: diffusive coeff = Gam    |
!     |_______________________________________________|

!     ----------------- 
      call grandientLSM(dudx,dudy,duds,uf,No_cp,nbe,sig) 
      call grandientLSM(dvdx,dvdy,dvds,vf,No_cp,nbe,sig) 
      call grandientLSM(dwdx,dwdy,dwds,wf,No_cp,nbe,sig) 
!     ----------------- 
      Cs  = 0.1d0
      do k=1,NZ
         do i=1,N_CELL  
            Sdudx = dudx(i,k)+sigmax(i,k)*duds(i,k)
            Sdvdx = dvdx(i,k)+sigmax(i,k)*dvds(i,k)
            Sdwdx = dwdx(i,k)+sigmax(i,k)*dwds(i,k)
            Sdudy = dudy(i,k)+sigmay(i,k)*duds(i,k)
            Sdvdy = dvdy(i,k)+sigmay(i,k)*dvds(i,k)
            Sdwdy = dwdy(i,k)+sigmay(i,k)*dwds(i,k)
            Sdudz = sigmaz(i,k)*duds(i,k)
            Sdvdz = sigmaz(i,k)*dvds(i,k)
            Sdwdz = sigmaz(i,k)*dwds(i,k)
!           ---------------
            DD  = ((0.5d0**2)*dsig(k)*Hpr(i))**(1.0d0/3.0d0)
!           ---------------
            nut = Hpr(i)*((Cs*DD)**2)
            Sij = (Sdudx+Sdudx)/2.0d0
            Gamux(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdudy+Sdvdx)/2.0d0    
            Gamuy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdudz+Sdwdx)/2.0d0
            Gamuz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
!           ----------
            Sij = (Sdvdx+Sdudy)/2.0d0
            Gamvx(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdvdy+Sdvdy)/2.0d0
            Gamvy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdvdz+Sdwdy)/2.0d0
            Gamvz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
!           ----------
            Sij = (Sdwdx+Sdudz)/2.0d0
            Gamwx(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdwdy+Sdvdz)/2.0d0
            Gamwy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            Sij = (Sdwdz+Sdwdz)/2.0d0
            Gamwz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  2.2)  Righ-hand side values                  |
!     |_______________________________________________|

!     ________________________________________________
!     Hydro case: p_hydro + gravity

      if (DynamicPressure.eq.0) then
         do k=1,NZ
            do i=1,N_CELL  
               rhsu(i,k) = (-gra*detadx(i))!*Hpr(i)
               rhsv(i,k) = (-gra*detady(i))!*Hpr(i) 
               rhsw(i,k) = 0.0d0          
            enddo
         enddo
!     ________________________________________________
!     Dynamic case: p_hydro + gravity + p_dynamic(old)
      else
!        -----------------
         call grandientLSM(dpdx,dpdy,dpds,pf,No_cp,nbe,sig) 
!        ----------------- 
         do k=1,NZ
            do i=1,N_CELL  
               Sdpdx = dpdx(i,k)+sigmax(i,k)*dpds(i,k)
               Sdpdy = dpdy(i,k)+sigmay(i,k)*dpds(i,k)
               Sdpdz = sigmaz(i,k)*dpds(i,k)
               rhsu(i,k) = (-gra*detadx(i)-1.0d0/rhof(i,k)*Sdpdx)*Hpr(i)
               rhsv(i,k) = (-gra*detady(i)-1.0d0/rhof(i,k)*Sdpdy)*Hpr(i)
               rhsw(i,k) = (-1.0d0/rhof(i,k)*Sdpdz)*Hpr(i)
            enddo
         enddo
      endif
!      _______________________________________________
!     |                                               |
!     |  2.3)  Vertical velocity: omegaf              |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL 
            omegaf(i,k) = sigmat(i,k)+uf(i,k)*sigmax(i,k) &
                                     +vf(i,k)*sigmay(i,k) &
                                     +wf(i,k)*sigmaz(i,k)
         enddo
      enddo

!      _______________________________________________
!     |                                               |
!     |  2.4) Old velocity components: Huf,Hvf,Hwf    |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL  
            Huf(i,k) = Hpr(i)*uf(i,k)
            Hvf(i,k) = Hpr(i)*vf(i,k)
            Hwf(i,k) = Hpr(i)*wf(i,k)
         enddo
      enddo
      do k=1,NZ-1
         do nv=1,N_VERT
            Hufv(nv,k) = Hprv(nv)*Hufv(nv,k)
            Hvfv(nv,k) = Hprv(nv)*Hvfv(nv,k)
            Hwfv(nv,k) = Hprv(nv)*Hwfv(nv,k)
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  2.5) Update velocity components: uf,vf,wf^(*)|
!     |_______________________________________________|

!     _________________________________________________
!     Update advection-diffusion equation

      tagBCu = 0
      tagBCv = 0
      tagBCw = 0
!     ----------------------------- 
      call FS_AdvDiff(ufstar,Hufv,                            & 
                      rhsu,Gamux,Gamuy,Gamuz,uf,vf,omegaf,pf, &  
                      Huf,xc,yc,sig,dsig,No_cp,nbe,           &
                      Hufv,xv,yv,sigv,dsigv,No_vp,nbev,       &
                      Hpr,h,eta,                              &
                      Hprv,hv,etav,                           &
                      tagBCu)
!     ----------------------------- 
      call FS_AdvDiff(vfstar,Hvfv,                            & 
                      rhsv,Gamvx,Gamvy,Gamvz,uf,vf,omegaf,pf, & 
                      Hvf,xc,yc,sig,dsig,No_cp,nbe,           &
                      Hvfv,xv,yv,sigv,dsigv,No_vp,nbev,       &
                      Hpr,h,eta,                              &
                      Hprv,hv,etav,                           &
                      tagBCv)
!     ----------------------------- 
      call FS_AdvDiff(wfstar,Hwfv,                            & 
                      rhsw,Gamwx,Gamwy,Gamwz,uf,vf,omegaf,pf, & 
                      Hwf,xc,yc,sig,dsig,No_cp,nbe,           &
                      Hwfv,xv,yv,sigv,dsigv,No_vp,nbev,       &
                      Hpr,h,eta,                              &
                      Hprv,hv,etav,                           &
                      tagBCw)

!     _________________________________________________
!     BC of w in the vertical direction
!     -----------------------------
!     Cell-center
      do i=1,N_CELL
         uB = 0.5d0*(ufstar(i,1)+ufstar(i,2))
         vB = 0.5d0*(vfstar(i,1)+vfstar(i,2))
         dwfdtB(i) = -uB*dhdx(i)-vB*dhdy(i)
         !wfstar(i,1) = dwfdtB(i)
!        -----------------
         uT = 0.5d0*(ufstar(i,NZ-1)+ufstar(i,NZ))
         vT = 0.5d0*(vfstar(i,NZ-1)+vfstar(i,NZ))
         dwfdtT(i) = Hpr(i)*detadt(i)+uT*detadx(i)+vT*detady(i)
         !wfstar(i,NZ) = dwfdtT(i)
      enddo
!     -----------------------------
!     Vertex
      call interpolation2D(dwfvdtB,xv,yv,No_vp,nbev,dwfdtB,xc,yc,No_cp,nbe)
      call interpolation2D(dwfvdtT,xv,yv,No_vp,nbev,dwfdtT,xc,yc,No_cp,nbe)
      do nv=1,N_VERT
         !Hwfv(nv,1)    = dwfvdtB(nv)
         !Hwfv(nv,NZ-1) = dwfvdtT(nv)
      enddo
!      _______________________________________________
!     |                                               |
!     |  2.6)        --- OPTIONAL CASES  ---          |
!     |_______________________________________________|

!     _________________________________________________
!     New values: Huf,Hvf,Hwf
      do k=1,NZ
         do i=1,N_CELL
            !Huf(i,k) = ufstar(i,k)
            !Hvf(i,k) = vfstar(i,k)
            !Hwf(i,k) = wfstar(i,k)
         enddo
      enddo
!     _________________________________________________
!     New values: omegaf
      do k=1,NZ
         do i=1,N_CELL
            !uf(i,k) = ufstar(i,k)/Hpr(i)
            !vf(i,k) = vfstar(i,k)/Hpr(i)
            !wf(i,k) = wfstar(i,k)/Hpr(i)
            !omegaf(i,k) = sigmat(i,k)+uf(i,k)*sigmax(i,k) &
            !                         +vf(i,k)*sigmay(i,k) &
            !                         +wf(i,k)*sigmaz(i,k)
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  2.7) --- OPTIONAL TEST: Exact velocity  ---  |
!     |_______________________________________________|

      UseThis = 0
      if (UseThis.eq.1) then
         do i=1,N_CELL
            x = xc(i)
            y = yc(i)
            etaE = FS_funeta(x,y,time)
            do k=1,NZ
               z = sig(k)*(etaE+h(i))-h(i)
               ufstar(i,k) = FS_funu(x,y,z,time)
               vfstar(i,k) = FS_funv(x,y,z,time)
               wfstar(i,k) = FS_funw(x,y,z,time)
            enddo
         enddo
      endif

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |   ----  3) ONLY HYDROSTATIC PRESSURE SOLUTION ----     |
!     |                    °  pfnp                             |
!     |                    °  ufnp                             |
!     |                    °  vfnp                             |    
!     |                    °  wfnp                             |        
!     |________________________________________________________|
!      ********************************************************

     IF (DynamicPressure.eq.0) THEN
         do i=1,N_CELL
            do k=1,NZ
               ufnp(i,k) = ufstar(i,k)/Hpr(i)
               vfnp(i,k) = vfstar(i,k)/Hpr(i)
               wfnp(i,k) = wfstar(i,k)/Hpr(i)
               pfnp(i,k) = pa + gra*rho_f*(1.d0-sig(k))*Hpr(i)
            enddo
         enddo
!        ----------------------------------------
!        Interpolation
         call interpolation3D(ufv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              ufnp,xc,yc,sig,dsig,No_cp,nbe)
         call interpolation3D(vfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              vfnp,xc,yc,sig,dsig,No_cp,nbe)
         call interpolation3D(wfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              wfnp,xc,yc,sig,dsig,No_cp,nbe)
         call interpolation3D(pfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              pfnp,xc,yc,sig,dsig,No_cp,nbe)

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |      ----  4) POISSON PRESSURE CALCULATION   ----      |
!     |               (only with dynamic pressure)             |        
!     |                    °  pf^(new)                         |
!     |                                                        |
!     |       ---- 5) VELOCITY CORRECTION STEP:  ----          |
!     |                    °  uf^(new)                         | 
!     |                    °  vf^(new)                         |   
!     |                    °  wf^(new)                         |
!     |________________________________________________________|
!      ********************************************************

      ELSE
!      _______________________________________________
!     |                                               |
!     |  4.1)  Coefficients: diffusive coeff = H*Gamp |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL  
            Gamp(i,k) = 1.0d0/rhof(i,k)
            Gamp(i,k) = Hpr(i)*Gamp(i,k)
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  4.2)  Pressure by projection method:         |
!     |        Dp = pfnp - pf                         |
!     |_______________________________________________|

      if (PressureMethod.eq.1) then
!        ______________________________________________
!        rhs
!        ----------------- 
         call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig) 
         call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig) 
         call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig) 
!        ----------------- 
         do k=1,NZ
            do i=1,N_CELL  
               Sdudx = dudx(i,k)+sigmax(i,k)*duds(i,k)
               Sdvdy = dvdy(i,k)+sigmay(i,k)*dvds(i,k)
               Sdwdz = sigmaz(i,k)*dwds(i,k)
               c1 = -ufstar(i,k)/Hpr(i)*dHprdx(i)
               c2 = -vfstar(i,k)/Hpr(i)*dHprdy(i)
               rhsp(i,k) = (Sdudx+Sdvdy+Sdwdz+c1+c2)/dt
            enddo
         enddo
!        __________________________________________________
!        Initial guess 
         do k=1,NZ
            do i=1,N_CELL  
               Dpfnp(i,k) = pfnp(i,k) - pf(i,k) 
            enddo
         enddo
!        __________________________________________________
!        Poisson solver  
!        ---------------
!        Variable Dp
         call FS_Poisson(Dpfnp,Dpfv,                  &
                         rhsp,Gamp,Gamp,Gamp,         & 
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,eta,                   &
                         Hprv,hv,etav,                &
                         rhof,dwfdtB,                 &
                         rhofv,dwfvdtB)
!        __________________________________________________
!        Pressure solution  
         do k=1,NZ
            do i=1,N_CELL  
               pfnp(i,k) = Dpfnp(i,k) + pf(i,k) 
            enddo
         enddo
         do k=1,NZ-1
            do nv=1,N_VERT  
               pfv(nv,k) = Dpfv(nv,k) + pfv(nv,k) 
            enddo
         enddo
!        __________________________________________________
!        Velocity correction
!        ----------------- 
         call grandientLSM(dDpdx,dDpdy,dDpds,Dpfnp,No_cp,nbe,sig) 
!        ----------------- 
         do k=1,NZ
            do i=1,N_CELL0
               SdDpdx = dDpdx(i,k)+sigmax(i,k)*dDpds(i,k)
               SdDpdy = dDpdy(i,k)+sigmay(i,k)*dDpds(i,k)
               SdDpdz = sigmaz(i,k)*dDpds(i,k)
               ufnp(i,k) = ufstar(i,k)-(dt/rhof(i,j))*SdDpdx
               vfnp(i,k) = vfstar(i,k)-(dt/rhof(i,j))*SdDpdy
               wfnp(i,k) = wfstar(i,k)-(dt/rhof(i,j))*SdDpdz 
            enddo
         enddo
      endif
!      _______________________________________________
!     |                                               |
!     |  4.3)  Pressure by direct method              |
!     |_______________________________________________|

      if (PressureMethod.eq.2) then
!        ______________________________________________
!        New velocity values: u* = u/H
         do k=1,NZ
            do i=1,N_CELL
               ufstar(i,k) = ufstar(i,k)/Hpr(i)
               vfstar(i,k) = vfstar(i,k)/Hpr(i)
               wfstar(i,k) = wfstar(i,k)/Hpr(i)
            enddo
         enddo
         do k=1,NZ-1
            do nv=1,N_VERT
               ufv(nv,k) = Hufv(nv,k)/Hprv(nv)
               vfv(nv,k) = Hvfv(nv,k)/Hprv(nv)
               wfv(nv,k) = Hwfv(nv,k)/Hprv(nv)
            enddo
         enddo
!        _______________________________________________
!        rhs
!        ----------------- 
         call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig) 
         call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig) 
         call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig) 
!        ----------------- 
         do k=1,NZ
            do i=1,N_CELL 
               Sdudx = dudx(i,k)+sigmax(i,k)*duds(i,k)
               Sdvdx = dvdx(i,k)+sigmax(i,k)*dvds(i,k)
               Sdwdx = dwdx(i,k)+sigmax(i,k)*dwds(i,k)
               Sdudy = dudy(i,k)+sigmay(i,k)*duds(i,k)
               Sdvdy = dvdy(i,k)+sigmay(i,k)*dvds(i,k)
               Sdwdy = dwdy(i,k)+sigmay(i,k)*dwds(i,k)
               Sdudz = sigmaz(i,k)*duds(i,k)
               Sdvdz = sigmaz(i,k)*dvds(i,k)
               Sdwdz = sigmaz(i,k)*dwds(i,k)

               rhsph = -gra*(d2etadx2(i)+d2etady2(i))       
               rhsp(i,k) = rhsph -2.0d0*(Sdudy*Sdvdx-Sdudx*Sdvdy &
                                        +Sdudz*Sdwdx-Sdudx*Sdwdz &
                                        +Sdvdz*Sdwdy-Sdvdy*Sdwdz)
               rhsp(i,k) = Hpr(i)*rhsp(i,k)
            enddo
         enddo
!        ______________________________________________
!        Initial guess
         do k=1,NZ
            do i=1,N_CELL  
               pfnp(i,k) = pf(i,k)
            enddo
         enddo
!        ______________________________________________
!        Pressure solver: pf^(n+1)
         call FS_Poisson(pfnp,pfv,                    &
                         rhsp,Gamp,Gamp,Gamp,         & 
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,eta,                   &
                         Hprv,hv,etav,                &
                         rhof,dwfdtB,                 &
                         rhofv,dwfvdtB)
!        ______________________________________________
!        Diffusion coefficients: Gam
         do k=1,NZ
            do i=1,N_CELL  
               Sdudx = dudx(i,k)+sigmax(i,k)*duds(i,k)
               Sdvdx = dvdx(i,k)+sigmax(i,k)*dvds(i,k)
               Sdwdx = dwdx(i,k)+sigmax(i,k)*dwds(i,k)
               Sdudy = dudy(i,k)+sigmay(i,k)*duds(i,k)
               Sdvdy = dvdy(i,k)+sigmay(i,k)*dvds(i,k)
               Sdwdy = dwdy(i,k)+sigmay(i,k)*dwds(i,k)
               Sdudz = sigmaz(i,k)*duds(i,k)
               Sdvdz = sigmaz(i,k)*dvds(i,k)
               Sdwdz = sigmaz(i,k)*dwds(i,k)
!              ---------------
               DD  = ((0.5d0**2)*dsig(k)*Hpr(i))**(1.0d0/3.0d0)
!              ---------------
               nut = Hpr(i)*((Cs*DD)**2)
               Sij = (Sdudx+Sdudx)/2.0d0
               Gamux(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdudy+Sdvdx)/2.0d0    
               Gamuy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdudz+Sdwdx)/2.0d0
               Gamuz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
!              ----------
               Sij = (Sdvdx+Sdudy)/2.0d0
               Gamvx(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdvdy+Sdvdy)/2.0d0
               Gamvy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdvdz+Sdwdy)/2.0d0
               Gamvz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
!               ----------
               Sij = (Sdwdx+Sdudz)/2.0d0
               Gamwx(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdwdy+Sdvdz)/2.0d0
               Gamwy(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdwdz+Sdwdz)/2.0d0
               Gamwz(i,k) = nut*dsqrt(2.0d0*Sij*Sij)
            enddo
         enddo
!        ______________________________________________
!        Right-hand side: rhs
!        -----------------
         call grandientLSM(dpdx,dpdy,dpds,pfnp,No_cp,nbe,sig) 
!        ----------------- 
         do k=1,NZ
            do i=1,N_CELL  
               Sdpdx = dpdx(i,k)+sigmax(i,k)*dpds(i,k)
               Sdpdy = dpdy(i,k)+sigmay(i,k)*dpds(i,k)
               Sdpdz = sigmaz(i,k)*dpds(i,k)
               rhsu(i,k) = (-gra*detadx(i)-1.0d0/rhof(i,k)*Sdpdx)*Hpr(i)
               rhsv(i,k) = (-gra*detady(i)-1.0d0/rhof(i,k)*Sdpdy)*Hpr(i)
               rhsw(i,k) = (-1.0d0/rhof(i,k)*Sdpdz)*Hpr(i)
            enddo
         enddo
!        ______________________________________________
!        Update velocity components
         tagBCu = 0
         tagBCv = 0
         tagBCw = 0
!        ----------------------------- 
         call FS_AdvDiff(ufnp,Hufv,                                & 
                         rhsu,Gamux,Gamuy,Gamuz,uf,vf,omegaf,pfnp, &  
                         Huf,xc,yc,sig,dsig,No_cp,nbe,             &
                         Hufv,xv,yv,sigv,dsigv,No_vp,nbev,         &
                         Hpr,h,eta,                                &
                         Hprv,hv,etav,                             &
                         tagBCu)
!        ----------------------------- 
         call FS_AdvDiff(vfnp,Hvfv,                                & 
                         rhsv,Gamvx,Gamvy,Gamvz,uf,vf,omegaf,pfnp, & 
                         Hvf,xc,yc,sig,dsig,No_cp,nbe,             &
                         Hvfv,xv,yv,sigv,dsigv,No_vp,nbev,         &
                         Hpr,h,eta,                                &
                         Hprv,hv,etav,                             &
                         tagBCv)
!        ----------------------------- 
         call FS_AdvDiff(wfnp,Hwfv,                                & 
                         rhsw,Gamwx,Gamwy,Gamwz,uf,vf,omegaf,pfnp, & 
                         Hwf,xc,yc,sig,dsig,No_cp,nbe,             &
                         Hwfv,xv,yv,sigv,dsigv,No_vp,nbev,         &
                         Hpr,h,eta,                                &
                         Hprv,hv,etav,                             &
                         tagBCw)
!        _________________________________________________
!        BC for w in the vertical direction
!        -----------------------------
!        Cell-center
         do i=1,N_CELL0
            uB = 0.5d0*(ufnp(i,1)+ufnp(i,2))
            vB = 0.5d0*(vfnp(i,1)+vfnp(i,2))
            dwfdtB(i) = -uB*dhdx(i)-vB*dhdy(i)
            !wfnp(i,1) = dwfdtB(i)
!           -----------------
            uT = 0.5d0*(ufnp(i,NZ-1)+ufnp(i,NZ))
            vT = 0.5d0*(vfnp(i,NZ-1)+vfnp(i,NZ))
            dwfdtT(i) = Hpr(i)*detadt(i)+uT*detadx(i)+vT*detady(i)
            !wfnp(i,NZ) = dwfdtT(i)
         enddo
!        -----------------------------
!        Vertex
         call interpolation2D(dwfvdtB,xv,yv,No_vp,nbev,dwfdtB,xc,yc,No_cp,nbe)
         call interpolation2D(dwfvdtT,xv,yv,No_vp,nbev,dwfdtT,xc,yc,No_cp,nbe)
         do nv=1,N_VERT
            !Hwfv(nv,1)    = dwfvdtB(nv)
            !Hwfv(nv,NZ-1) = dwfvdtT(nv)
         enddo
      endif
!     _________________________________________________
!     New values: ufnp,vfnp,wfnp
      do k=1,NZ
         do i=1,N_CELL
            ufnp(i,k) = ufnp(i,k)/Hpr(i)
            vfnp(i,k) = vfnp(i,k)/Hpr(i)
            wfnp(i,k) = wfnp(i,k)/Hpr(i)
         enddo
      enddo
      do k=1,NZ-1
         do nv=1,N_VERT
            ufv(nv,k) = Hufv(nv,k)/Hprv(nv)
            vfv(nv,k) = Hvfv(nv,k)/Hprv(nv)
            wfv(nv,k) = Hwfv(nv,k)/Hprv(nv)
         enddo
      enddo

      ENDIF !<--- dynamic pressure

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
                 dudx,dudy,duds,        &
                 dvdx,dvdy,dvds,        &
                 dwdx,dwdy,dwds,        &
                 dpdx,dpdy,dpds,        &
                 dDpdx,dDpdy,dDpds,     &
                 Gamux,Gamuy,Gamuz,     &
                 Gamvx,Gamvy,Gamvz,     &
                 Gamwx,Gamwy,Gamwz,     &
                 Gamp,                  &
                 rhsu,rhsv,rhsw,rhsp,   &
                 ph,omegaf,             &
                 Huf,Hvf,Hwf,           &
                 Hufv,Hvfv,Hwfv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TEST Free Surface'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          END OF FREE SURFACE TEST WITH NAVIER-STOKES EQUATION       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
