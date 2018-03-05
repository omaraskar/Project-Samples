!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     VERTEX BOUNDARY CONDITION                       !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                            phi2D,xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition.   !
!    We have three choices:                                           !
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
!  | <-- phiv2D  |   (N_VERT)  | Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!  | * i,j,nv    |   Loop counters                                 |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
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

      real*8, dimension(:)  :: phiv2D(N_VERT)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
!     ---------------------------------
      real*8, dimension(:)  :: phi2D(N_CELL)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:) :: dphidx(N_CELL)
      real*8 ,dimension(:) :: dphidy(N_CELL)
      real*8 :: sumfx,sumfy,deter
!     ---------------------------------
      real*8 :: nnx,nny
      real*8 :: x1,y1,x2,y2,z1,z2,z3
!     -----------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: x,y,fB,dfBdn
      real*8 :: Neumanndfdn2D
      real*8 :: funSolExam2D
!     ---------------------------------
      integer :: ic1,ic2,jv1,jv2,jv3,jj,jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: VertexBC2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Vertex Dirichlet Boundary Condition                !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.1) THEN
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
               x = xv(nv)
               y = yv(nv)
               phiv2D(nv) = funSolExam2D(x,y)
            endif
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!                  Vertex Neumann BC approximation                    !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.2) THEN
!      ________________________________________________________
!     |                                                        |
!     |           Approximation of the gradient 2D             |
!     |________________________________________________________|

      do i=1,N_CELL0	
	 sumfx = 0.0d0
	 sumfy = 0.0d0
	 do j=1,3
	    jc = No_cp(i,j)                
	    sumfx = sumfx + dxCC(i,j)*(phi2D(jc)-phi2D(i))
	    sumfy = sumfy + dyCC(i,j)*(phi2D(jc)-phi2D(i))
	 enddo
	 deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
         dphidx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	 dphidy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                        Neumann BC                      |  
!     |________________________________________________________|

      do nv=1,N_VERT
         if (nbev(nv).ne.0) then
!           ---------------------------
!           Exact Neumann BC
            x = xv(nv)
            y = yv(nv)
            nnx = normxv(nv)
            nny = normyv(nv)
            dfBdn = -Neumanndfdn2D(x,y,nnx,nny) 
!           ---------------------------
!           Function approximations
            ic1 = icxn1(nv)
            ic2 = icxn2(nv)
            f1 = phi2D(ic1) + dphidx(ic1)*(xn1(nv)-xc(ic1)) &
                            + dphidy(ic1)*(yn1(nv)-yc(ic1))   
            f2 = phi2D(ic2) + dphidx(ic2)*(xn2(nv)-xc(ic2)) &
                            + dphidy(ic2)*(yn2(nv)-yc(ic2))
!           ---------------------------
!           Approx at the vertex point
            h1 = 1.0d0*dn(nv)
            h2 = 2.0d0*dn(nv)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv2D(nv) = f0 
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
           print*, '      <----   End subroutine: VertexBC2D'
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
