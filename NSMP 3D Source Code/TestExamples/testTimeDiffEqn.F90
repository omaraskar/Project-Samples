!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   TEST ADVECTION-DIFFUSION EQUATION                 !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeDiffEqn(phiA,phivA,                       &  
                                 phi,xc,yc,sig,dsig,No_cp,nbe,     &
                                 phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 nRK)              
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the Diffusion equation in each time step.   !
!                      d(phi)/dt  = DIFF + rhs                        ! 
!                                                                     !
!    We have four time discretization methods to solve the problem.   !
!    (chosen at file cppdefs.h).                                      ! 
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name      |    Size     | Description                     |  !  
!  |_______________|_____________|_________________________________|  !
!  | <-- phiA      |(N_CELL,NZ)  | Approximate solution cell-center|  !
!  | <-- phiE      |(N_CELL,NZ)  | Exact solution cell-center      |  !  
!  | <-- Error     |(N_CELL,NZ)  | Error cell-center               |  !    
!  | <-- phivA     |(N_VERT,NZ-1)| Approximate solution vertex     |  !
!  | <-- phivE     |(N_VERT,NZ-1)| Exact solution vertex           |  ! 
!  | <-- Errorv    |(N_VERT,NZ-1)| Error vertex                    |  !  
!  |_______________|_____________|_________________________________|  !
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

      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phivA(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      integer :: nRK
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:) :: phiA2D(N_CELL)
      real*8,dimension(:) :: phivA2D(N_VERT)
      real*8,dimension(:) :: phi2D(N_CELL)
      real*8,dimension(:) :: phiv2D(N_VERT)
!     ---------------------------------------- 
      real*8,dimension(:) :: uu2D(N_CELL)
      real*8,dimension(:) :: vv2D(N_CELL)
      real*8,dimension(:) :: Gamx2D(N_CELL)
      real*8,dimension(:) :: Gamy2D(N_CELL)
      real*8,dimension(:) :: rhs2D(N_CELL)
!     -------------------------------------   
      real*8,dimension(:,:) :: uu(N_CELL,NZ) 
      real*8,dimension(:,:) :: vv(N_CELL,NZ) 
      real*8,dimension(:,:) :: ww(N_CELL,NZ) 
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
      real*8,dimension(:,:) :: rhs(N_CELL,NZ) 
!     -------------------------------------
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
!     ----------------------------------------
      real*8 :: TimeExample2D,TimeExample3D
      real*8 :: Peakvalue,PeakvalueEx,PeakError
      real*8 :: x,y,z
      integer:: DisplayThis
!     --------------------------------------
      integer:: tagBC
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Test DiffEqn'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                       Diffusion approximation  2D                   !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.2) THEN
!        ________________________________________________________
!        2D assignation 

         do i=1,N_CELL
            phi2D(i) = phi(i,2)
         enddo
         do nv=1,N_VERT
            phiv2D(nv) = phiv(nv,2)
         enddo
!         ________________________________________________________
!        |                                                        |
!        |               Components of the problem                |
!        |________________________________________________________|

         do i=1,N_CELL
!           -----------------------------
!           Velocity profile
            uu2D(i) = 0.0d0
            vv2D(i) = 0.0d0   
!           -----------------------------
!           Diffusive coefficient
            Gamx2D(i) = 5.0d-2 
            Gamy2D(i) = 5.0d-2 
!           -----------------------------
!           Right-hand side
            rhs2D(i) = 0.0d0 
         enddo
!         ________________________________________________________
!        |                                                        |
!        |                        Solution                        |
!        |________________________________________________________|

         tagBC=1
         call AdvDiffEqn2D(phiA2D,phivA2D,                &
                           uu2D,vv2D,Gamx2D,Gamy2D,rhs2D, & 
                           phi2D,xc,yc,No_cp,nbe,         &
                           phiv2D,xv,yv,No_vp,nbev,       &
                           tagBC)
!         ________________________________________________________
!        |                                                        |
!        |                     Final solution                     |
!        |________________________________________________________|

         do k=1,NZ 
            do i=1,N_CELL0
               phiA(i,k) = phiA2D(i)
            enddo
         enddo
         do k=1,NZ-1 
            do nv=1,N_VERT
               phivA(nv,k) = phivA2D(nv)
            enddo
         enddo

!*********************************************************************!
!                                                                     !
!                       Diffusion approximation  3D                   !
!                                                                     !
!*********************************************************************!

      ELSEIF (TestDimension.eq.3) THEN 
!         ________________________________________________________
!        |                                                        |
!        |               Components of the problem                |
!        |________________________________________________________|

         do k=1,NZ
            do i=1,N_CELL
!              -----------------------------
!              Velocity profile
               uu(i,k) = 0.0d0
               vv(i,k) = 0.0d0 
               ww(i,k) = 0.0d0  
!              -----------------------------
!              Diffusive coefficient
               Gamx(i,k) = 5.0d-2 
               Gamy(i,k) = 5.0d-2 
               Gamz(i,k) = 5.0d-2 
!              -----------------------------
!              Right-hand side
               rhs(i,k) = 0.0d0                
            enddo
         enddo
!         ________________________________________________________
!        |                                                        |
!        |                       Solution                         |
!        |________________________________________________________|

         tagBC=1
         call AdvDiffEqn3D(phiA,phivA,                       &
                           uu,vv,ww,Gamx,Gamy,Gamz,rhs,      &  
                           phi,xc,yc,sig,dsig,No_cp,nbe,     &
                           phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           tagBC)
      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Test DiffEqn'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF testAdvection                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
