!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               INITIALIZATION OF THE FREE-SURFACE PROBLEM            !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_testTimeInitial(ufn,vfn,wfn,pfn,             &
                                    ufv,vfv,wfv,pfv,             &
                                    xc,yc,sig,dsig,No_cp,nbe,    &
                                    xv,yv,sigv,dsigv,No_vp,nbev, & 
                                    Hpr,h,etan,                  &
                                    Hprv,hv,etav,                &
                                    xct,yct,zct,                 & 
                                    xvt,yvt,zvt,                 &
                                    No_sp)             
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the initial condition of the time test   !  
!    problems with free surface.                                      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |           Name       |    Size     | Description              |  !  
!  |______________________|_____________|__________________________|  !
!  | <--- ufn,vfn,wfn,pfn |(N_CELL,NZ)  | Initial solution center  |  !
!  | <--- ufv,vfv,wfv,pfv |(N_VERT,NZ-1)| Initial solution vertex  |  !
!  |______________________|_____________|__________________________|  !
!  | <--- Hpr,etan        |(N_CELL)     | Initial solution center  |  !
!  | <--- Hprv,etavv      |(N_VERT)     | Initial solution vertex  |  !
!  |______________________|_____________|__________________________|  !
!  | <--- xct,yct,zct     |(N_CELL,NZ)  | Coordinates center       |  !
!  | <--- xvt,yvt,zvt     |(N_VERT,NZ-1)| Coordinates vertex       |  !
!  |______________________|_____________|__________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc |(N_CELL)    | Coordinates of the cell centers     |  !
!  | ---> sig   |(NZ)        | Sigma value at the cell centers     |  !
!  | ---> dsig  |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | ---> No_cp |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | ---> nbe   |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | ---> sigv  |(NZ-1)      | sigma of the vertex points          |  !
!  | ---> dsigv |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
!  | ---> No_vp |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | ---> nbev  |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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

      real*8,dimension(:,:) :: ufn(N_CELL,NZ)
      real*8,dimension(:,:) :: vfn(N_CELL,NZ)
      real*8,dimension(:,:) :: wfn(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
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
      real*8,dimension(:,:) :: xct(N_CELL,NZ)
      real*8,dimension(:,:) :: yct(N_CELL,NZ)
      real*8,dimension(:,:) :: zct(N_CELL,NZ)
      real*8,dimension(:,:) :: xvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: yvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: zvt(N_VERT,NZ-1)
!     --------------------------------------
      integer,dimension(:)  :: No_sp(N_SPmax)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: FS_funu,FS_funv,FS_funw,FS_funeta
      real*8 :: MaxErrorA,sumErrorA,x,y,z

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: FS_TimeInitial'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Initial time of the Example                     !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                       Free surface                     |
!     |________________________________________________________|

!     -----------------------------------
!     Cell-center   
      do i=1,N_CELL
         x = xc(i)
         y = yc(i)
         etan(i) = FS_funeta(x,y,time)
         Hpr(i)  = etan(i) + h(i)
      enddo
!     -----------------------------------
!     Vertex  
      do nv=1,N_VERT
         x = xv(nv)
         y = yv(nv)
         etav(nv) = FS_funeta(x,y,time) 
         Hprv(nv) = etav(nv) + hv(nv)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Velocity and pressure                  |
!     |________________________________________________________|

!     -----------------------------------
!     Cell-center  
      do i=1,N_CELL
         do k=1,NZ
            x = xc(i)
            y = yc(i)
            z = sig(k)*Hpr(i)-h(i)
            ufn(i,k) = FS_funu(x,y,z,time)
            vfn(i,k) = FS_funv(x,y,z,time)
            wfn(i,k) = FS_funw(x,y,z,time)
            pfn(i,k) = 0.0d0 
         enddo
      enddo
!     -----------------------------------
!     Vertex 
      do nv=1,N_VERT
         do k=1,NZ-1
            x = xv(nv)
            y = yv(nv)
            z = sig(k)*Hprv(nv)-hv(nv)                                                  
    	    ufv(nv,k) = FS_funu(x,y,z,time) 
            vfv(nv,k) = FS_funv(x,y,z,time) 
            wfv(nv,k) = FS_funw(x,y,z,time) 
            pfv(nv,k) = 0.0d0 
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Coordinates (xt,yt,zt)                 |
!     |________________________________________________________|

!     -----------------------------------
!     Cell-center       
      do i=1,N_CELL
         do k=1,NZ
            xct(i,k) = xc(i)
            yct(i,k) = yc(i)
            zct(i,k) = sig(k)*Hpr(i)-h(i)
         enddo
      enddo  
!     -----------------------------------
!     Vertex 
      do nv=1,N_VERT
         do k=1,NZ-1
            xvt(nv,k) = xv(nv)
            yvt(nv,k) = yv(nv)
            zvt(nv,k) = sigv(k)*Hprv(nv)-hv(nv)
         enddo
      enddo 
!      ________________________________________________________
!     |                                                        |
!     |                 Save initial results                   |
!     |________________________________________________________|

!     ______________________________________________________
!     Write initial time error file 
      maxErrorA = 0.0d0
      sumErrorA = 0.0d0
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(9100,*) time,maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
      IF (rang_topo.eq.0) THEN
         write(9100,*) time,maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA,&
                            maxErrorA,sumErrorA
      ENDIF
#     endif
!     =============== END ================    
!     ====================================

!     ______________________________________________________
!     Write initial time sample solution 

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         nv = No_sp(1)
         write(9101,*) time,xv(nv),yv(nv),etav(nv),&
                            ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1)
         nv = No_sp(2)
         write(9102,*) time,xv(nv),yv(nv),etav(nv),&
                            ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1)
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
      IF (rang_topo.eq.0) THEN
         nv = No_sp(1)
         write(9101,*) time,xv(nv),yv(nv),etav(nv),&
                            ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1)
         nv = No_sp(2)
         write(9102,*) time,xv(nv),yv(nv),etav(nv),&
                            ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1)
      ENDIF
#     endif
!     =============== END ================    
!     ====================================

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: FS_TimeInitial'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                     END OF testAdvection Initial                    !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
