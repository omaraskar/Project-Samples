!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!              BOUNDARY CONDITIONS OF THE POISSON PROBLEM             !
!                              May 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCpoisson(phi,phiv,                        &        
                           xc,yc,sig,dsig,No_cp,nbe,        &
                           xv,yv,sigv,dsigv,No_vp,nbev)      

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine updates the boundary condition values of the     !
!    poisson variable for the new time step simulation.               ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--  phi   |(N_CELL,NZ)  | function phi at the cell center    |  !
!  | <--  phiv  |(N_VERT,NZ-1)| function phi at the vertex         |  !
!  |____________|_____________|____________________________________|  !
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
!  | <--funOldv |(N_VERT,NZ) | Cell-vertex solution at t(n)        |  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !  
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbe    |(N_CELL0)   | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !

!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- N_CELL0 | Number of inside cells                          |  !
!  |     NZ      | Points in the sigma direction                   |  !
!  |     time    | time                                            |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  | fB          | function at the edge point of a cell boundary   |  ! 
!  | dfBdn       | normal derivative of fB                         |  !
!  | jj,jv1,jv2  | index for cell neighbor and vertices            |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!   <-->  Input variables                                             !
!   ----  Parameters                                                  !
!         Common variables used                                       !
!    *    Common variables modified                                   !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ) 
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      real*8 :: fB,dfBdn
      integer:: jj,jv1,jv2,jv3,j1,j2,elem
!     ----------------------------------------
      real*8,dimension(:,:) :: pv(N_VERT,NZ-1) 
      real*8 :: x,y,z,uu,vv,ww
!     ----------------------------------------
      real*8 :: funEta,funh
!     ----------------------------------------
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: funExample1u,funExample1v,funExample1w
      real*8 :: funExample1Dp,funExample1p

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: BC poisson'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                   Sigma transform: Depth: H = h + eta               !
!                                                                     !
!*********************************************************************!
!     ________________________________________
!     Cell-center 
      do i=1,N_CELL 
         etan(i) = funEta(xc(i),yc(i))
         h(i)    = funh(xc(i),yc(i)) 
         Hpr(i)  = h(i) + etan(i)  
      enddo
!     ________________________________________
!     Vertex
      do nv=1,N_VERT 
         etav(nv) = funEta(xv(nv),yv(nv)) 
         hv(nv)   = funh(xv(nv),yv(nv))
         Hprv(nv) = hv(nv) + etav(nv)
      enddo
!     ________________________________________
!     Box case
      if (ChooseDomBox.eq.1) then
         do i=1,N_CELL 
            etan(i) = 0.0d0
            h(i)    = 0.0d0
            Hpr(i)  = 1.0d0  
         enddo
         do nv=1,N_VERT 
            etav(nv) = 0.0d0
            hv(nv)   = 0.0d0
            Hprv(nv) = 1.0d0
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                 Boundary conditions  cell-centers                   !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                       Horizontal                       |
!     |________________________________________________________|


      DO i=1,N_CELL0
!        ______________________________________________________
!        Wall
         IF (nbe(i).eq.1) THEN
 	    do j=1,3
!              -------------------------------
!              Numbering index cell-cente
	       nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================
!                 -------------------------------
!                 Vertex of the boundary (not used now)
		  jj=j+1
		  if (jj.gt.3) jj=jj-3
		  jv1 = No_vp(i,j)
		  jv2 = No_vp(i,jj)
!                 --------------------------------
!                 Outside point value
                  if (TestNSProblem.eq.1) then
                     do k=1,NZ   
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExample1Dp(x,y,z,time)                                                   
    	                phi(nc,k) = 2.0d0*fB-phi(i,k) 
                     enddo
                  endif
	       endif
            enddo
!        ______________________________________________________
!        Discharge boundary
         ELSEIF (nbe(i).eq.2) THEN
 	    do j=1,3
!              -------------------------------
!              Numbering index cell-cente
	       nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================
!                 -------------------------------
!                 Vertex of the boundary (not used now)
		  jj=j+1
		  if (jj.gt.3) jj=jj-3
		  jv1 = No_vp(i,j)
		  jv2 = No_vp(i,jj)
!                 --------------------------------
!                 Outside point value
                  if (TestNSProblem.eq.1) then
                     do k=1,NZ   
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExample1Dp(x,y,z,time)                                                   
    	                phi(nc,k) = 2.0d0*dlCE(i,j)*fB + phi(i,k)  
                     enddo
                  endif
	       endif
            enddo
!        ______________________________________________________
!        Water level
         ELSEIF (nbe(i).eq.3) THEN
 	    do j=1,3
!              -------------------------------
!              Numbering index cell-cente
	       nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================    
!                 ====================================
!                 -------------------------------
!                 Vertex of the boundary (not used now)
		  jj=j+1
		  if (jj.gt.3) jj=jj-3
		  jv1 = No_vp(i,j)
		  jv2 = No_vp(i,jj)
!                 --------------------------------
!                 Outside point value
                  if (TestNSProblem.eq.1) then
                     do k=1,NZ   
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExample1Dp(x,y,z,time)                                                   
    	                phi(nc,k) = 2.0d0*dlCE(i,j)*fB+phi(i,k) 
                     enddo
                  endif
	       endif
            enddo
         ENDIF
      ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                         Vertical                       |
!     |________________________________________________________|

      IF (ChooseBoundary.eq.1) THEN
         if (TestDimension.eq.3) then
!           _______________________________
!           Example: N-S problem 
            if (TestNSProblem.eq.1) then
               do i=1,N_CELL
                  x = xc(i)
                  y = yc(i)
!                 ------------------------ 
!                 Bottom
                  z  = sig(1)*Hpr(i)-h(i)
	          fB = funExample1Dp(x,y,z) 
                  phi(i,1) = 2.0d0*fB-phi(i,2) 
!                 ------------------------ 
!                 Top
                  z = sig(NZ)*Hpr(i)-h(i)
	          fB = funExample1Dp(x,y,z) 
                  phi(i,NZ) = 2.0d0*fB-phi(i,NZ-1)

               enddo
            endif
         endif
      ELSEIF (ChooseBoundary.eq.2) THEN
         if (TestDimension.eq.3) then
!           _______________________________
!           Example: N-S problem 
            if (TestNSProblem.eq.1) then
               do i=1,N_CELL
                  x = xc(i)
                  y = yc(i)
!                 ------------------------ 
!                 Bottom
                  z  = sigv(1)*Hpr(i)-h(i)
	          fB = funExample1Dp(x,y,z) 
                  phi(i,1) = dsig(1)*fB+phi(i,2) 
!                 ------------------------ 
!                 Top
                  z = sigv(NZ-1)*Hpr(i)-h(i)
	          fB = funExample1Dp(x,y,z) 
                  phi(i,NZ) = dsig(NZ)*fB+phi(i,NZ-1)
               enddo
            endif
         endif
      ENDIF

!*********************************************************************!
!                                                                     !
!                   Boundary conditions at the vertices               !
!  (If we use Neumann then we should obtain them from cell-centers)   !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.1) THEN
         if (TestNSProblem.eq.1) then
!           _______________________________
!           Dirichlet
            if (nbev(nv).eq.1) then
               do nv=1,N_VERT
                  x = xv(nv)
                  y = yv(nv)
                  do k=1,NZ-1                                                   
                     z = sigv(k)
                     phiv(nv,k) = funExample1p(x,y,z,time)  
                  enddo
               enddo
            endif
         endif
      ENDIF

!       call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                            phi,xc,yc,sig,dsig,No_cp)


!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: BC poisson'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	          END                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
