!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      Finalization one time step                     !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE boundaryCondition(phi,                             &        
                                   xc,yc,sig,dsig,No_cp,nbe,        &
                                   xv,yv,sigv,dsigv,No_vp,nbev)      

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine updates the boundary condition values of the     !
!    main variables for the new time step simulation.                 ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> phi   |(N_CELL,NZ) | function phi                        |  !
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

      real*8 :: fB,dfBdn
      integer:: jj,jv1,jv2,elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: Boundary conditions '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         Boundary conditions                         !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                       Horizontal                       |
!     |________________________________________________________|

      DO i=1,N_CELL0
!        _______________________________
!       |                               |
!       |             Wall              |
!       |_______________________________|

         IF (nbe(i).eq.1) THEN
 	    do j=1,3
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
!                 -------------------------------
!                 Value at the edge point 
		  fB = 0.0d0 
!                 -------------------------------
!                 Outside point value    
                  do k=1,NZ                                                   
    	             phi(nc,k) = 2.0d0*fB - phi(i,k) 
                  enddo
	       endif
            enddo
!        _______________________________
!       |                               |
!       |      Discharge boundary       |
!       |_______________________________|

         ELSEIF (nbe(i).eq.2) THEN
 	    do j=1,3
!              -------------------------------
!              Numbering index cell-center  
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
!                 Numbering index vertex 
		  jj=j+1
		  if (jj.gt.3) jj=jj-3
		  jv1 = No_vp(i,j)
		  jv2 = No_vp(i,jj)
!                 -------------------------------
!                 Value at the edge point 
                  fB = 0.0d0                          !888888888888888888888888
		  dfBdn = 0.0d0
!                 -------------------------------
!                 Outside point value
                  do k=1,NZ                                                   
    	             phi(nc,k) = dlCC(i,j)*dfBdn + phi(i,k)                      
    	             phi(nc,k) = 2.0d0*fB - phi(i,k)  !888888888888888888888888    
                  enddo			
	       endif
            enddo
!        _______________________________
!       |                               |
!       |          Water level          |
!       |_______________________________|

         ELSEIF (nbe(i).eq.3) THEN
 	    do j=1,3
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
!                 -------------------------------
!                 Value at the edge point 
                  fB = 0.0d0                         !888888888888888888888888
		  dfBdn = 0.0d0
!                 -------------------------------
!                 Outside point value   
                  do k=1,NZ                                                   
    	             phi(nc,k) = dlCC(i,j)*dfBdn + phi(i,k) 
    	             phi(nc,k) = 2.0d0*fB - phi(i,k)  !888888888888888888888888     
                  enddo					
	       endif
            enddo
         ENDIF
      ENDDO

!      ________________________________________________________
!     |                                                        |
!     |                         Vertical                       |
!     |________________________________________________________|

      if (TestDimension.eq.3) then
         do i=1,N_CELL
	    phi(i,1)  = 0.0d0   ! Bottom    
	    phi(i,NZ) = 0.0d0   ! Top 
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: Boundary conditions'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	          END                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
