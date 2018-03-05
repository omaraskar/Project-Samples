!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!        2D BOUNDARY CONDITION CELL-CENTERS OF THE FREE SURFACE       !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_BCeta(Hnew,h,xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates Dirichlet cell-center boundary           !
!    condition for the free surface.                                  !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- Hnew   |   (N_CELL)  | Function at the vertices          |  !
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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

      real*8, dimension(:)  :: Hnew(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
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
      real*8 :: FS_funeta
      integer :: elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: FS_BCeta'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Dirichlet Boundary Condition                    !
!                                                                     !
!*********************************************************************!

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
                  fB = FS_funeta(x,y,time) + h(i)
                  Hnew(nc) = 2.0d0*fB-Hnew(i)
               endif 
	    enddo
         endif
      enddo
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: FS_BCeta'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END of Free Surface BC                      !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
