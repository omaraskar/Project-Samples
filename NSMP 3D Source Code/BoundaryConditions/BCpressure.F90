!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                 PRESSURE BOUNDARY CONDITION CELL-CENTERS            !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCprescenter3D(phi,xc,yc,sig,dsig,No_cp,nbe)
  
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition for the pressure. We have two choices:                 !
!    BoundaryType = 1:  All boundaries of the problem are Dirichlet   !
!    BoundaryType = 2:  All boundaries of the proble are Neumann      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- phi     |(N_CELL,NZ)| Function at the cell-center         |  !   
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB,dfBdn
      real*8 :: nnx,nny,nnz
      real*8 :: funExamNSp,NeumanndpdnNS
!     --------------------------------------
      integer, parameter :: BoundaryType = 1 
!     --------------------------------------
      integer :: elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCprescenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Dirichlet Boundary Condition                    !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Dirichlet                       |
!     |________________________________________________________|

      IF (BoundaryType.eq.1) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = sig(k)+0.5d0*dsig(k)
            fB = funExamNSp(x,y,z,time) 
            phi(i,k)  = 2.0d0*fB-phi(i,k+1)
!           ______                   
!           Top
            k = NZ
            z = sig(k)-0.5d0*dsig(k)
            fB = funExamNSp(x,y,z,time)
            phi(i,k)= 2.0d0*fB-phi(i,k-1)
!           --------------------------------------------------
!           HORIZONTAL
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
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        fB = funExamNSp(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                        Neumann                         |
!     |________________________________________________________|

      ELSEIF (BoundaryType.eq.1) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = sig(k)+0.5d0*dsig(k)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0
            dfBdn = 0.0d0!NeumanndpdnNS(x,y,z,nnx,nny,nnz)
            phi(i,k) = dsig(k)*dfBdn + phi(i,k+1) 
!           ______                   
!           Top
            k = NZ
            z = sig(k)-0.5d0*dsig(k)
            nnx = 0.0d0
            nny = 0.0d0
            nnz = 1.0d0
	    dfBdn = 0.0d0!NeumanndpdnNS(x,y,z,nnx,nny,nnz)
            phi(i,k) = dsig(k-1)*dfBdn + phi(i,k-1)
!           --------------------------------------------------
!           HORIZONTAL
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
                     nnx = normxc(i,j)
                     nny = normyc(i,j)
                     nnz = 0.0d0
                     do k=1,NZ 
                        x = xe(i,j)
                        y = ye(i,j)
                        z = sig(k)
                        dfBdn = 0.0d0!NeumanndpdnNS(x,y,z,nnx,nny,nnz)
                        phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn + phi(i,k)
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
      ENDIF
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCprescenter3D'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF CELL-CENTER BOUNDARY CONDITION            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
