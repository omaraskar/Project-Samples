!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     VERTEX BOUNDARY CONDITION                       !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                            phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
                             
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
!  | <-- phiv    |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!  | * i,j,nv,k  |   Loop counters                                 |  !    
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
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  | --> Hprv    |(N_VERT)   | Total depth = h + etan  vertices    |  !
!  | --> hv      |(N_VERT)   | still depth             vertices    |  !
!  | --> etav    |(N_VERT)   | free surface            vertices    |  !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     ---------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:,:):: dphidx(N_CELL,NZ)
      real*8 ,dimension(:,:):: dphidy(N_CELL,NZ)
      real*8 :: sumfx,sumfy,deter
!     ---------------------------------
      real*8 :: nnx,nny,nnz
      integer :: ic1,ic2
!     -----------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: x,y,z,fB,dfBdn
      real*8 :: funSolExam3D,Neumanndfdn3D
      real*8 :: funExamNSp,NeumanndpdnNS
!     ---------------------------------
      integer :: jv1,jv2,jv3,jj,jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: VertexBC3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                        Navier-Stokes Equation                       !
!                                                                     !
!*********************************************************************!

      IF (RUNTestNSEqn.eq.1) THEN

!      ________________________________________________________
!     |                                                        |
!     |               Dirichlet Boundary Condition             |
!     |________________________________________________________|

      if (ChooseBoundary.eq.1) then
!        ________________________________________________________
!        Horizontal
         do k=1,NZ-1 
            do nv=1,N_VERT
               if (nbev(nv).ne.0) then
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)
                  phiv(nv,k) = funExamNSp(x,y,z,time)
               endif
            enddo
         enddo
!        ________________________________________________________
!        Vertical
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)
            phiv(nv,k) = funExamNSp(x,y,z,time)
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)  
            phiv(nv,k) = funExamNSp(x,y,z,time)
         enddo
      endif

!      ________________________________________________________
!     |                                                        |
!     |                 Neumann BC approximation               |
!     |________________________________________________________|

      if (ChooseBoundary.eq.2) then
!        ________________________________________________________
!        Horizontal

!        --------------------------------------------------------  
!        Approximation of the gradient 2D      
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
!        ---------------------------------------------------------  
!        Neumann BC approximation  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
!              ---------------------------
!              Exact Neumann BC
               x = xv(nv)
               y = yv(nv)
               nnx = normxv(nv)
               nny = normyv(nv)
               do k=1,NZ-1
                  z = sigv(k)
                  nnz =  0.0d0
                  dfBdn = -NeumanndpdnNS(x,y,z,nnx,nny,nnz,time)
!                 ---------------------------
!                 Function approximations
                  ic1 = icxn1(nv)
                  ic2 = icxn2(nv)
                  f1 = phi(ic1,k) + dphidx(ic1,k)*(xn1(nv)-xc(ic1)) &
                                  + dphidy(ic1,k)*(yn1(nv)-yc(ic1))   
                  f2 = phi(ic2,k) + dphidx(ic2,k)*(xn2(nv)-xc(ic2)) &
                                  + dphidy(ic2,k)*(yn2(nv)-yc(ic2))
!                 ---------------------------
!                 Approx at the vertex point
                  h1 = 1.0d0*dn(nv)
                  h2 = 2.0d0*dn(nv)
                  deno = (h1*h2*h2-h2*h1*h1)
                  f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
                  phiv(nv,k) = f0 
               enddo
            endif           
         enddo
!        ________________________________________________________
!        Vertical

         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
!           _________                    
!           Bottom
            k  = 1
            z  = sigv(k)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0
            dfBdn = -NeumanndpdnNS(x,y,z,nnx,nny,nnz,time) 
            f1 = phiv(nv,k+1) 
            f2 = phiv(nv,k+2) 
            h1 = sigv(k+1)-sigv(k)
            h2 = sigv(k+2)-sigv(k)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv(nv,k) = f0 
!           _________                   
!           Top
            k = NZ-1
            z = sigv(k)  
            nnx =  0.0d0
            nny =  0.0d0
            nnz =  1.0d0
            dfBdn = NeumanndpdnNS(x,y,z,nnx,nny,nnz,time)
            f1 = phiv(nv,k-1) 
            f2 = phiv(nv,k-2) 
            h1 = sigv(k-1)-sigv(k)
            h2 = sigv(k-2)-sigv(k)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv(nv,k)= f0  
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                            Other problems                           !
!                                                                     !
!*********************************************************************!
      ELSE
!      ________________________________________________________
!     |                                                        |
!     |               Dirichlet Boundary Condition             |
!     |________________________________________________________|

      if (ChooseBoundary.eq.1) then
!        ________________________________________________________
!        Horizontal
         do k=1,NZ-1 
            do nv=1,N_VERT
               if (nbev(nv).ne.0) then
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  phiv(nv,k) = funSolExam3D(x,y,z) 
               endif
            enddo
         enddo
!        ________________________________________________________
!        Vertical
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k=1
            z  = sigv(k)*Hprv(nv)-hv(nv)
            phiv(nv,k) = funSolExam3D(x,y,z) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)*Hprv(nv)-hv(nv)
            phiv(nv,k) = funSolExam3D(x,y,z)
         enddo
      endif

!      ________________________________________________________
!     |                                                        |
!     |                 Neumann BC approximation               |
!     |________________________________________________________|

      if (ChooseBoundary.eq.2) then
!        ________________________________________________________
!        Horizontal

!        --------------------------------------------------------  
!        Approximation of the gradient 2D      
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
!        ---------------------------------------------------------  
!        Neumann BC approximation  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
!              ---------------------------
!              Exact Neumann BC
               x = xv(nv)
               y = yv(nv)
               nnx = normxv(nv)
               nny = normyv(nv)
               do k=1,NZ-1
                  z = sigv(k)
                  nnz =  0.0d0
                  dfBdn = -Neumanndfdn3D(x,y,z,nnx,nny,nnz)
!                 ---------------------------
!                 Function approximations
                  ic1 = icxn1(nv)
                  ic2 = icxn2(nv)
                  f1 = phi(ic1,k) + dphidx(ic1,k)*(xn1(nv)-xc(ic1)) &
                                  + dphidy(ic1,k)*(yn1(nv)-yc(ic1))   
                  f2 = phi(ic2,k) + dphidx(ic2,k)*(xn2(nv)-xc(ic2)) &
                                  + dphidy(ic2,k)*(yn2(nv)-yc(ic2))
!                 ---------------------------
!                 Approx at the vertex point
                  h1 = 1.0d0*dn(nv)
                  h2 = 2.0d0*dn(nv)
                  deno = (h1*h2*h2-h2*h1*h1)
                  f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
                  phiv(nv,k) = f0 
               enddo
            endif           
         enddo
!        ________________________________________________________
!        Vertical

         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
!           _________                    
!           Bottom
            k  = 1
            z  = sigv(k)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0
            dfBdn = -Neumanndfdn3D(x,y,z,nnx,nny,nnz) 
            f1 = phiv(nv,k+1) 
            f2 = phiv(nv,k+2) 
            h1 = sigv(k+1)-sigv(k)
            h2 = sigv(k+2)-sigv(k)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv(nv,k) = f0 
!           _________                   
!           Top
            k = NZ-1
            z = sigv(k)  
            nnx =  0.0d0
            nny =  0.0d0
            nnz =  1.0d0
            dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
            f1 = phiv(nv,k-1) 
            f2 = phiv(nv,k-2) 
            h1 = sigv(k-1)-sigv(k)
            h2 = sigv(k-2)-sigv(k)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv(nv,k)= f0  
         enddo
      endif
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
      
      ENDIF
 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: VertexBC3D'
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
