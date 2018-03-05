!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  NEUMANN VERTEX BOUNDARY CONDITION                  !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE NeumannVertexBC(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 phi,xc,yc,sig,dsig,No_cp,nbe)

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
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8,dimension(:,:) :: sigv(NZ-1)
      real*8,dimension(:,:) :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
!     ---------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8,dimension(:,:) :: sig(NZ)
      real*8,dimension(:,:) :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:,:):: dphidx(N_CELL,NZ)
      real*8 ,dimension(:,:):: dphidy(N_CELL,NZ)
      real*8 :: sumfx,sumfy,deter
!     ---------------------------------
      real*8 :: dn,nnx,nny,nnz
      real*8 :: xn1,yn1,xn2,yn2
      real*8 :: x1,y1,x2,y2,z1,z2,z3
!     -----------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: x,y,z,fB,dfBdn
      real*8 :: Neumanndfdn2D,Neumanndfdn3D
      real*8 :: funSolExam2D,funSolExam3D
!     ---------------------------------
      integer :: jv1,jv2,jv3,jj,jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: VertexBC'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Vertex Dirichlet Boundary Condition                !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.1) THEN
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
            if (nbev(nv).ne.0) then
               phiv2D(nv) = funSolExam2D(x,y)
            endif
         enddo
      ENDIF

      IF (ChooseBoundary.eq.1) THEN
         do i=1,N_CELL0
            do k=2,NZ-1 
               x = xc(i)
               y = yc(i)
               z = sig(k)
               phi(i,k) = funSolExam3D(x,y,z) 
            enddo
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!                  Vertex Neumann BC approximation                    !
!                                                                     !
!*********************************************************************!

      IF (ChooseBoundary.eq.1) THEN
!      ________________________________________________________
!     |                                                        |
!     |           Approximation of the gradient 2D             |
!     |________________________________________________________|

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

!      ________________________________________________________
!     |                                                        |
!     |                        Neumann BC                      |  
!     |________________________________________________________|


      DO i=1,N_CELL0
!        --------------------------------------------------
!        Length dn in the normal direction
         dn = 0.45d0*min(dlVV(i,1),dlVV(i,2),dlVV(i,3))
 	 do j=1,3
	    nv=No_vp(i,j)
	    if (nbev(nv).ne.0) then  
!              _____________________________________________
!              Unit normal n in the out direction
               nnx = normxv(nv)
               nny = normyv(nv)
               !print*,i,'nv=',nv,'nnx,nny=',nnx,nny
!              _____________________________________________
!              Points in the normal direction
               xn1 = xv(nv)+dn*(-nnx)
               yn1 = yv(nv)+dn*(-nny)
               xn2 = xv(nv)+2.0d0*dn*(-nnx)
               yn2 = yv(nv)+2.0d0*dn*(-nny)
               !print*,'            nv=',nv,'(xn2,yn2)=',xn2,yn2
!              _____________________________________________
!              Is the point inside the triangle?
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)

               x1 = dxVV(i,1)
               y1 = dyVV(i,1)
               x2 = xn2-xv(jv1)
               y2 = yn2-yv(jv1)
               z1 = x1*y2-y1*x2

               x1 = dxVV(i,2)
               y1 = dyVV(i,2)
               x2 = xn2-xv(jv2)
               y2 = yn2-yv(jv2)
               z2 = x1*y2-y1*x2

               x1 = dxVV(i,3)
               y1 = dyVV(i,3)
               x2 = xn2-xv(jv3)
               y2 = yn2-yv(jv3)
               z1 = x1*y2-y1*x2
               !print*,'            nv=',nv,'z1,z2,z3',z1,z2,z3
!              _____________________________________________
!              If YES, approximate the vertex point
               if ((z1.ge.0).and.(z2.ge.0).and.(z3.ge.0)) then
                  !print*,'            nv=',nv,'YES'
                  do k=1,NZ-1
!                    ---------------------------
!                    Exact Neumann BC
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)
                     nnz =  0.0d0
                     !dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz) 
                     dfBdn = -Neumanndfdn2D(x,y,nnx,nny) 
!                    ---------------------------
!                    Function approximations
                     f1 = phi(i,k) + dphidx(i,k)*(xn1-xc(i)) &
                                   + dphidy(i,k)*(yn1-yc(i))   
                     f2 = phi(i,k) + dphidx(i,k)*(xn2-xc(i)) &
                                   + dphidy(i,k)*(yn2-yc(i))
!                    ---------------------------
!                    Approx at the vertex point
                     h1 = 1.0d0*dn
                     h2 = 2.0d0*dn
                     deno = (h1*h2*h2-h2*h1*h1)
                     f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)&
                           /(h2*h2-h1*h1)
                     phiv(nv,k) = f0
                  enddo 
                  !print*,'            nv=',nv,'fn1,fn2=',f1,f2
                  !print*,'------------------------'
               else
                  !print*,'            nv=',nv,'NO'
                  !print*,'------------------------'
               endif           
	    endif
         enddo
      ENDDO

      ENDIF
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: VertexBC'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                      	   END OF INTERPOLATION                       !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
