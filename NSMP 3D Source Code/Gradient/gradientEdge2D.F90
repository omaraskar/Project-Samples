!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             CALCULATION OF THE GRADIENT AT EACH FACE 2D             !
!                              March 2013                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientEdge2D(ax,ay,az,bx,by,bz,fx,fy,fz,  &
                                 fun,xc,yc,sig,dsig,No_cp,    &
                                 funv,xv,yv,sigv,dsigv,No_vp)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program approximates the gradient of a function fun at      !
!    the 3 edges of a triangle taking an approximation of the         !
!    integral over a 2D region surronding each edge (Appendix B).     !
!    The program output are the corresponding constants next to       !
!    the average values:                                              !
!                                                                     !
!                   d(fun)/dx= ax*fun(i)+bxfun(j)+fx                  !
!                   d(fun)/dy= ay*fun(i)+byfun(j)+fy                  !
!                                                                     !
!    where fx,fy correspond to the values of the fun at the vert-     !
!    tices of the triangle.                                           ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |    Size      |     Description                 |  !  
!  |______________|______________|_________________________________|  !
!  | <-- ax,bx,fx |(N_CELL0,NZ,5)|  d(fun)/dx coefficients         |  !
!  | <-- ay,by,fy |(N_CELL0,NZ,5)|  d(fun)/dy coefficients         |  !   
!  |______________|______________|_________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> fun     |(N_CELL,NZ)| Function "fun" at the center element|  !
!  | --> xc,yc   |(N_CELL)   | Cell-center corrdinates             |  !
!  | --> sig     |(NZ)       | Vertical location of the cell center|  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!  | --> funv    |(N_VERT,NZ-1)| Function "fun" at the vertices    |  !
!  | --> xv,yv   |(N_VERT)   | Vertex corrdinates                  |  !
!  | --> sigv    |(NZ-1)     | Vertical location of the vertices   |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k+1)               |  !
!  | --> No_vp   |(N_VERT,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  |--- N_CELL   |  Total number of the cells                      |  !
!  |--- N_CELL0  |  Number of the cell centers inside the domain   |  !
!  |    N_VERT   |  Number of the vertices in the horizontal domain|  !
!  |    NZ       |  Number of vertical points                      |  ! 
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  !   
!  | *nv,nc,i,j,k|  Loop counters: vertices,cells, other           |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !  
!  |  sumfx      | real variable to add gradient component x       |  !
!  |  sumfy      | real variable to add gradient component y       |  !
!  |  jv1,jv2    | Number index of the vertices                    |  !
!  |  s          | Loop counter of the number of the region edges  |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |   Keys and common parameters                           |
!     |________________________________________________________|

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
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of variables                            |
!     |________________________________________________________|

      real*8,dimension(:,:,:):: ax(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: ay(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: az(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: bx(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: by(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: bz(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: fx(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: fy(N_CELL0,NZ,5)
      real*8,dimension(:,:,:):: fz(N_CELL0,NZ,5)

      real*8,dimension(:,:) :: fun(N_CELL,NZ)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)

      real*8,dimension(:,:) :: funv(N_VERT,NZ-1)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1) 
      integer,dimension(:,:):: No_vp(N_CELL0,3)

!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________| 

      real*8,dimension(:) :: nxDLoArea(1:4)
      real*8,dimension(:) :: nyDLoArea(1:4)
      real*8 :: funvp
      real*8 :: sumax,sumay,sumfx
      real*8 :: sumbx,sumby,sumfy
      real*8 :: sumax1,sumbx1,sumfx1
      real*8 :: AreaReg
      integer:: jc,s,kkk
      integer:: jj,jv1,jv2,jv3
      integer,parameter :: ChooseOptionEdge = 2

      real*8 :: temfv1,temfv2
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientEdge2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF (ChooseOptionEdge.eq.1) THEN

!*********************************************************************!
!                                                                     !
!               Calculation of the gradient: Option 1                 !
!                                                                     !
!*********************************************************************!

      DO k=1,NZ
         do i=1,N_CELL0
            do j=1,3
!              _______________________________
!             |                               |
!             |             Index             |
!             |_______________________________|

!              --------------------------------
!              Cell-center index
	       jc = No_cp(i,j)
!              --------------------------------
!              Vertex index               
	       jj=j+1
	       if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              --------------------------------
!              Volume of the region               
               AreaReg = (areaCell(i)+areaCell(jc))/3.0d0 
!              --------------------------------
!              nx*Length
	       nxDLoArea(1) = (yc(jc)-yv(jv1))/AreaReg
	       nxDLoArea(2) = (yv(jv2)-yc(jc))/AreaReg
	       nxDLoArea(3) = (yc(i)-yv(jv2))/AreaReg
	       nxDLoArea(4) = (yv(jv1)-yc(i))/AreaReg
!              --------------------------------
!              ny*length
	       nyDLoArea(1) = -(xc(jc)-xv(jv1))/AreaReg
	       nyDLoArea(2) = -(xv(jv2)-xc(jc))/AreaReg
	       nyDLoArea(3) = -(xc(i)-xv(jv2))/AreaReg
	       nyDLoArea(4) = -(xv(jv1)-xc(i))/AreaReg

!              _______________________________
!             |                               |
!             |  Edges of the region: s=1,2   |
!             |_______________________________|

	       sumbx = 0.0d0
	       sumby = 0.0d0
               sumfx = 0.0d0
	       sumfy = 0.0d0
               do s=1,2
!                 -------------------------------
!                 Function at the vertices
                  if (s.eq.1) then 
                     funvp = funv(jv1,k)
                  elseif (s.eq.2) then 
                     funvp = funv(jv2,k)
                  endif
!                 -------------------------------
!                 Constant cell neighbor j 
                  sumbx=sumbx+0.5d0*nxDLoArea(s)
                  sumby=sumby+0.5d0*nyDLoArea(s)
!                 -------------------------------
!                 Known vertex values              
	          sumfx=sumfx+0.5d0*nxDLoArea(s)*funvp
	          sumfy=sumfy+0.5d0*nyDLoArea(s)*funvp
               enddo !s

!              _______________________________
!             |                               |
!             |  Edges of the region: s=3,4   |
!             |_______________________________|

               sumax = 0.0d0
	       sumay = 0.0d0
               do s=3,4
!                 -------------------------------
!                 Function at the vertices
                  if (s.eq.3) then 
                     funvp = funv(jv2,k)
                  elseif (s.eq.4) then 
                     funvp = funv(jv1,k)
                  endif
!                 -------------------------------
!                 Constant cell-center i 
                  sumax=sumax+0.5d0*nxDLoArea(s)
                  sumay=sumay+0.5d0*nyDLoArea(s)
!                 -------------------------------
!                 Known vertex values              
	          sumfx=sumfx+0.5d0*nxDLoArea(s)*funvp
	          sumfy=sumfy+0.5d0*nyDLoArea(s)*funvp
               enddo !s
!              _______________________________
!             |                               |
!             |         Contributions         |
!             |_______________________________|

!              ----------------------------
!              Constant of cell center i   
               ax(i,k,j) = sumax
	       ay(i,k,j) = sumay
!              ----------------------------
!              Constant of neighbors j    
               bx(i,k,j) = sumbx
	       by(i,k,j) = sumby
!              ----------------------------
!              Constant of known vertex    
               fx(i,k,j) = sumfx
	       fy(i,k,j) = sumfy
	    enddo !j
         enddo
      ENDDO
 

      ELSEIF (ChooseOptionEdge.eq.2) THEN

!*********************************************************************!
!                                                                     !
!               Calculation of the gradient: Option 2                 !
!                                                                     !
!*********************************************************************!

         DO k=1,NZ
            do i=1,N_CELL0
               do j=1,3
!                 --------------------------------
!                 Cell-center index
	          jc = No_cp(i,j)
!                 --------------------------------
!                 Vertex index               
	          jj=j+1
	          if (jj.gt.3) jj=jj-3
	          jv1 = No_vp(i,j)
	          jv2 = No_vp(i,jj)
!                 --------------------------------
!                 Volume of the region               
                  AreaReg = (areaCell(i)+areaCell(jc))/3.0d0
!                 --------------------------------
!                 Constants x
                  sumax = -0.5d0*(yv(jv2)-yv(jv1))/AreaReg
                  sumbx =  0.5d0*(yv(jv2)-yv(jv1))/AreaReg
                  sumfx = -0.5d0*(yc(jc)-yc(i))/AreaReg
!                 --------------------------------
!                 Constants y
                  sumay =  0.5d0*(xv(jv2)-xv(jv1))/AreaReg
                  sumby = -0.5d0*(xv(jv2)-xv(jv1))/AreaReg
                  sumfy =  0.5d0*(xc(jc)-xc(i))/AreaReg
!                 --------------------------------
!                 Constant of cell center i   
                  ax(i,k,j) = sumax
	          ay(i,k,j) = sumay
!                 --------------------------------
!                 Constant of neighbors j    
                  bx(i,k,j) = sumbx
	          by(i,k,j) = sumby
!                 --------------------------------
!                 Known vertex values       
                  fx(i,k,j) = sumfx*(funv(jv2,k)-funv(jv1,k))
	          fy(i,k,j) = sumfy*(funv(jv2,k)-funv(jv1,k)) 
	       enddo 
            enddo
         ENDDO

      ENDIF

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientEdge2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                        END GRADIENT EDGE 2D                         !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
