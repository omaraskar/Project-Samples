!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   PRESSURE BOUNDARY CONDITIONS FOR FREE SURFACE PROBLEMS (center)   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_BCpc(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan,&
                         rhof,dwfdtB)
  
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the boundary condition for the pressure  !
!    in the case of free surface problems at the cell-centers.        !
!                                                                     !      
!    In the vertical walls:                                           !
!                              dpdn = 0               (Neumann)       !
!    At the free surface:                                             !
!                                 p = 0               (Dirichlet)     !
!    At the bottom boundary:                                          !
!                            dpdsig = -rho*H*dwdt     (Neumann)       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--> phi    |(N_CELL,NZ)| Function at the cell-center         |  !   
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
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !  
!  | --> rhof    |(N_CELL,NZ)| Density                             |  !
!  | --> dwdtB   |(N_CELL)   | d(w)/dt at the bottom boundary      |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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
!     --------------------------------------
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8,dimension(:,:) :: rhof(N_CELL,NZ)
      real*8,dimension(:,:) :: dwfdtB(N_CELL)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: fB,dfBdn,dwdt
      integer :: elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: FS_BCpc'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                            Initilaization                           !
!                                                                     !
!*********************************************************************!

!     ________________________________________________________
!     Vertical
      do i=1,N_CELL0
!        ______                    
!        Bottom
	 dfBdn = -rhof(i,1)*Hpr(i)*dwfdtB(i)
         phi(i,1) = dsig(1)*dfBdn + phi(i,2) 
!        ______                   
!        Top
         fB  = 0.0d0
         phi(i,NZ) = 2.0d0*fB-phi(i,NZ-1)
      enddo
!     ________________________________________________________
!     Horizontal
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
                  do k=1,NZ 
                     dfBdn = 0.0d0
                     phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn + phi(i,k)
                  enddo
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
           print*, '      <----   End subroutine: FS_BCpc'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   PRESSURE BOUNDARY CONDITIONS FOR FREE SURFACE PROBLEMS (vertex)   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_BCpv(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                         phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan,&
                         rhofv,dwfvdtB)
                             
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the boundary condition for the pressure  !
!    in the case of free surface problems at the vertex points.       !
!                                                                     !      
!    In the vertical walls:                                           !
!                              dpdn = 0               (Neumann)       !
!    At the free surface:                                             !
!                                 p = 0               (Dirichlet)     !
!    At the bottom boundary:                                          !
!                            dpdsig = -rho*H*dwdt     (Neumann)       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <-> phiv    |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
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
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!  | --> Hprv    |(N_VERT)   | Total depth = h + etan  vertices    |  !
!  | --> hv      |(N_VERT)   | still depth             vertices    |  !
!  | --> etav    |(N_VERT)   | free surface            vertices    |  !
!  |_____________|___________|_____________________________________|  !
!  | --> rhofv   |(N_VERT,NZ)| Density at the vertices             |  !
!  | --> dwdtB   |(N_VERT)   | d(w)/dt at the bottom boundary      |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8,dimension(:,:) :: rhofv(N_VERT,NZ-1)
      real*8, dimension(:)  :: dwfvdtB(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:):: dphidx(N_CELL,NZ)
      real*8,dimension(:,:):: dphidy(N_CELL,NZ)
      real*8 :: sumfx,sumfy,deter
!     --------------------------------------
      integer:: ic1,ic2
!     --------------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: x,y,z,fB,dfBdn
!     --------------------------------------
      integer:: jv1,jv2,jv3,jj,jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: FS_BCpv'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                              Initialization                         !
!                                                                     !
!*********************************************************************!

!     ________________________________________________________
!     Vertical
      do nv=1,N_VERT
!        _________                 
!        Bottom
         k  = 1
         dfBdn = -rhofv(nv,k)*Hprv(nv)*dwfvdtB(nv) 
         f1 = phiv(nv,k+1) 
         f2 = phiv(nv,k+2) 
         h1 = sigv(k+1)-sigv(k)
         h2 = sigv(k+2)-sigv(k)
         deno = (h1*h2*h2-h2*h1*h1)
         f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
         phiv(nv,k) = f0 
!        _________                
!        Top
         phiv(nv,NZ-1) = 0.0d0
      enddo

!     ________________________________________________________
!     Horizontal

!     -------------------------------------------------  
!     Approximation of the gradient 2D      
      do k=1,NZ
         do i=1,N_CELL0	
	    sumfx = 0.0d0
	    sumfy = 0.0d0
	    do j=1,3
	       jc = No_cp(i,j)                
	       sumfx = sumfx + dxCC(i,j)*(phi(jc,k)-phi(i,k))
	       sumfy = sumfy + dyCC(i,j)*(phi(jc,k)-phi(i,k))
	    enddo
	    deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
            dphidx(i,k)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	    dphidy(i,k)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
         enddo
      enddo
!     -------------------------------------------------  
!     Neumann BC approximation  
      do nv=1,N_VERT
         if (nbev(nv).ne.0) then
            do k=1,NZ-1
               dfBdn = 0.0d0
!              -------------------------
!              Function approximations
               ic1 = icxn1(nv)
               ic2 = icxn2(nv)
               f1 = phi(ic1,k) + dphidx(ic1,k)*(xn1(nv)-xc(ic1)) &
                               + dphidy(ic1,k)*(yn1(nv)-yc(ic1))   
               f2 = phi(ic2,k) + dphidx(ic2,k)*(xn2(nv)-xc(ic2)) &
                               + dphidy(ic2,k)*(yn2(nv)-yc(ic2))
!              --------------------------
!              Approx at the vertex point
               h1 = 1.0d0*dn(nv)
               h2 = 2.0d0*dn(nv)
               deno = (h1*h2*h2-h2*h1*h1)
               f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
               phiv(nv,k) = f0 
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
           print*, '      <----   End subroutine: FS_BCpv'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF PRESSURE BOUNDARY CONDITION               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
