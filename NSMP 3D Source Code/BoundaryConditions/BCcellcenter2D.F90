!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  3D BOUNDARY CONDITION CELL-CENTERS                 !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition. We have three choices:                                !
!                                                                     !      
!    ChooseBoundary = 0:  Combination of type of boundary. It is      !
!                         determintated by nbe & nbev.                !
!    ChooseBoundary = 1:  All boundaries of the problem are Dirichlet !
!    ChooseBoundary = 2:  All boundaries of the proble are Neumann    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <-- phi2D   |   (N_CELL)  | Function at the vertices          |  !
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Local variables:                                                 !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | xee,yee     |  Perpendicular intersection of each edge        |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
!   ---  Parameters                                                   !
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

      real*8, dimension(:)  :: phi2D(N_CELL)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,fB,dfBdn
      real*8 :: nnx,nny
      real*8 :: Neumanndfdn2D
      real*8 :: funSolExam2D
      integer :: elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCcellcenter2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Dirichlet Boundary Condition                    !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.1) THEN
         do i=1,N_CELL0
	    if (nbe(i).ne.0) then	
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
                     x = xe(i,j)
                     y = ye(i,j)
                     fB = funSolExam2D(x,y)   
                     phi2D(nc) = 2.0d0*fB-phi2D(i)
                  endif 
	       enddo
            endif
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!                       Neumann BC approximation                      !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.2) THEN
         do i=1,N_CELL0
	    if (nbe(i).ne.0) then	
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
                      x = xe(i,j)
                      y = ye(i,j)
                      nnx = normxc(i,j)
                      nny = normyc(i,j)
                      dfBdn = Neumanndfdn2D(x,y,nnx,nny)
                      phi2D(nc) = 2.0d0*dlCE(i,j)*dfBdn+phi2D(i)
                  endif 
	        enddo
            endif
         enddo
      ENDIF
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCcellcenter2D'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF VERTEX BOUNDARY CONDITION                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
