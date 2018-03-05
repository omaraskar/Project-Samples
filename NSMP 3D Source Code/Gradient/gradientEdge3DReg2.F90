!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               CALCULATION OF THE GRADIENT AT EACH FACE              !
!                              March 2013                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientEdge3DReg2(ax,ay,az,bx,by,bz,fx,fy,fz,  &
                                     fun,xc,yc,sig,dsig,No_cp,    &
                                     funv,xv,yv,sigv,dsigv,No_vp)         
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program approximates the gradient of a function fun at      !
!    the 5 faces of the prism (horizontal: j=1,2,3, top: j=4 and      !
!    bottom: j=5) taking approximations of the integral over a        !
!    3D region of type 2 surronding each face (Appendix B).           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |    Size      |     Description                 |  !  
!  |______________|______________|_________________________________|  !
!  | <-- ax,bx,fx |(N_CELL0,NZ,5)| d(fun)/dx= ax*fun(i)+bxfun(j)+fx|  !
!  | <-- ay,by,fy |(N_CELL0,NZ,5)| d(fun)/dy= ay*fun(i)+byfun(j)+fy|  ! 
!  | <-- az,bz,fz |(N_CELL0,NZ,5)| d(fun)/dz= az*fun(i)+bzfun(j)+fz|  !   
!  |______________|______________|_________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> fun     |(N_CELL,NZ)| Function "fun" at the center centers|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  |_____________|___________|_____________________________________|  !
!  | --> funv    |(N_VERT,NZ)| Function "fun" at the vertices      |  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |--- N_VERT  |  Number of the vertices in the horizontal domain |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |____________|__________________________________________________|  !
!  |  dlVV      |(N_CELL0,3)| length of the triangle side vert-vert|  !
!  |  dhCE      |(N_CELL0,3)| perpendiculat distance center-edge   |  !
!  |____________|___________|______________________________________|  !
!                                                                     !
!    Common integers:                                                 !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |       Name        |               Description                 |  !  
!  |___________________|___________________________________________|  !  
!  | funv1T,funv1B     | function at the triangle vertex point 1   |  !
!  | funv2T,funv2B     | function at the triangle vertex point 2   |  !
!  | funv3T,funv3B     | function at the triangle vertex point 3   |  !
!  | funvT,funvB       | function at the triangle vertex top & bott|  !
!  | sumax,sumbx,sumfx | real variable to add gradient component x |  !
!  | sumay,sumby,sumfy | real variable to add gradient component y |  !
!  | sumaz,sumbz,sumfz | real variable to add gradient component z |  !
!  | VoluReg           | Volumeof the Region                       |  !
!  | nxAreaoVol        | = (normalx)*(AreaFace)/VolumeRegion       |  !
!  | nyAreaoVol        | = (normaly)*(AreaFace)/VolumeRegion       |  !
!  | nzAreaoVol        | = (normalz)*(AreaFace)/VolumeRegion       |  !
!  | a1,a2,a3          | First displacement vector X1=(a1,a2,a3)   |  !
!  | b1,b2,b3          | Second displacement vector X2=(b1,b2,b3)  |  !
!  | jv1,jv2,jv3       | Numbering index of the vertices           |  !
!  | jc                | Numbering index of the cell neighbor      |  !
!  | jj,ss             | Index related to the vertex               |  !
!  | s                 | Loop counter of the faces of the region   |  !
!  |___________________|___________________________________________|  !
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

      real*8 ,dimension(:,:,:) :: ax(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: ay(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: az(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: bx(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: by(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: bz(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: fx(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: fy(N_CELL0,NZ,5)
      real*8 ,dimension(:,:,:) :: fz(N_CELL0,NZ,5)

      real*8 ,dimension(:,:):: fun(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)

      real*8 ,dimension(:,:):: funv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1) 
      integer,dimension(:,:):: No_vp(N_CELL0,3)

!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________| 

      real*8 :: funvT,funv1T,funv2T,funv3T
      real*8 :: funvB,funv1B,funv2B,funv3B
      real*8 :: sumfx,sumfy,sumfz
      real*8 :: sumax,sumay,sumaz
      real*8 :: sumbx,sumby,sumbz
      real*8 :: VoluReg
      real*8 :: nxAreaoVol,nyAreaoVol,nzAreaoVol
      real*8 :: a1,a2,a3
      real*8 :: b1,b2,b3
      integer:: jc,s,ss
      integer:: jj,jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientEdgeReg2'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Calculation of the gradient                     !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1
         do i=1,N_CELL0

!           ____________________________________________________
!          |                                                    |
!          |            Horizontal neighbors:  j=1:3            |
!          |____________________________________________________|

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
               sumax = 0.0d0
	       sumay = 0.0d0
	       sumaz = 0.0d0
	       sumbx = 0.0d0
	       sumby = 0.0d0
	       sumbz = 0.0d0
               sumfx = 0.0d0
	       sumfy = 0.0d0
	       sumfz = 0.0d0
!              _______________________________
!             |                               |
!             |      Volume of the region     |
!             |_______________________________|

               VoluReg = dlVV(i,j)*dsigv(k-1)*(dhCE(i,j))/3.0d0 

!              _______________________________
!             |                               |
!             |       Region faces: 1,2       |
!             |_______________________________|

               do s=1,2
                  if (s.eq.1) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv1,k-1)
                     funvT = funv(jv1,k)
!                    ----------------------------
!                    normal*Area
                     nxAreaoVol = 0.5d0*(yc(jc)-yv(jv1))*dsigv(k-1)/VoluReg
                     nyAreaoVol =-0.5d0*(xc(jc)-xv(jv1))*dsigv(k-1)/VoluReg
                  elseif (s.eq.2) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv2,k-1)
                     funvT = funv(jv2,k)
!                    ----------------------------
!                    normal*Area
	             nxAreaoVol = 0.5d0*(yv(jv2)-yc(jc))*dsigv(k-1)/VoluReg
	             nyAreaoVol =-0.5d0*(xv(jv2)-xc(jc))*dsigv(k-1)/VoluReg
                  endif
!                 -------------------------------
!                 Constant cell neighbor j 
                  sumbx = sumbx + (1.0d0/3.0d0)*nxAreaoVol
                  sumby = sumby + (1.0d0/3.0d0)*nyAreaoVol
!                 -------------------------------
!                 Known vertex values              
	          sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*(funvB+funvT)
	          sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*(funvB+funvT)
               enddo !s
!              _______________________________
!             |                               |
!             |       Region faces: 3,4       |
!             |_______________________________|

               do s=3,4
                  if (s.eq.3) then 
!                    --------------------------------------------
!                    Function at the points
                     funvB = funv(jv2,k-1)
                     funvT = funv(jv2,k)
!                    ----------------------------
!                    normal*Area
	             nxAreaoVol =  0.5d0*(yc(i)-yv(jv2))*dsigv(k-1)/VoluReg
	             nyAreaoVol = -0.5d0*(xc(i)-xv(jv2))*dsigv(k-1)/VoluReg
                  elseif (s.eq.4) then 
!                    --------------------------------------------
!                    Function at the points
                     funvB = funv(jv1,k-1)
                     funvT = funv(jv1,k)
!                    ----------------------------
!                    normal*Area
	             nxAreaoVol =  0.5d0*(yv(jv1)-yc(i))*dsigv(k-1)/VoluReg
	             nyAreaoVol = -0.5d0*(xv(jv1)-xc(i))*dsigv(k-1)/VoluReg
                  endif
!                 -------------------------------
!                 Cell-center i 
                  sumax = sumax + (1.0d0/3.0d0)*nxAreaoVol
                  sumay = sumay + (1.0d0/3.0d0)*nyAreaoVol
!                 -------------------------------
!                 Known vertex values              
	          sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*(funvB+funvT)
	          sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*(funvB+funvT)
               enddo

!              _______________________________
!             |                               |
!             |       Region faces: 5,6       |
!             |_______________________________|

!              ----------------------------
!              Function at the points
               funv1T = funv(jv1,k)
               funv2T = funv(jv2,k)
               funv1B = funv(jv1,k-1)
               funv2B = funv(jv2,k-1)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(jc)
               a2 = yv(jv2)-yc(jc)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(jc)
               b2 = yv(jv1)-yc(jc)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
               sumbx = sumbx + (2.0d0/3.0d0)*nxAreaoVol
               sumby = sumby + (2.0d0/3.0d0)*nyAreaoVol
!              -------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               ((funv1T+funv2T)+(funv1B+funv2B))
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               ((funv1T+funv2T)+(funv1B+funv2B))
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               ((funv1T+funv2T)-(funv1B+funv2B))
!              _______________________________
!             |                               |
!             |       Region faces: 7,8       |
!             |_______________________________|

!              ----------------------------
!              Function at the points
               funv1T = funv(jv1,k)
               funv2T = funv(jv2,k)
               funv1B = funv(jv1,k-1)
               funv2B = funv(jv2,k-1)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
               sumax = sumax + (2.0d0/3.0d0)*nxAreaoVol
               sumay = sumay + (2.0d0/3.0d0)*nyAreaoVol
!              -------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               ((funv1T+funv2T)+(funv1B+funv2B))
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               ((funv1T+funv2T)+(funv1B+funv2B))
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               ((funv1T+funv2T)-(funv1B+funv2B))  
!              _______________________________
!             |                               |
!             |         Contributions         |
!             |_______________________________|

!              ----------------------------
!              Cell center i   
               ax(i,k,j) = sumax
	       ay(i,k,j) = sumay
               az(i,k,j) = sumaz
!              ----------------------------
!              Cell neighbor j    
               bx(i,k,j) = sumbx
	       by(i,k,j) = sumby
               bz(i,k,j) = sumbz
!              ----------------------------
!              Known vertex values      
               fx(i,k,j) = sumfx
	       fy(i,k,j) = sumfy
               fz(i,k,j) = sumfz

	    enddo !j

!           ____________________________________________________
!          |                                                    |
!          |                Vertical TOP neighbor               |
!          |____________________________________________________|

            sumax = 0.0d0
	    sumay = 0.0d0
	    sumaz = 0.0d0
            sumbx = 0.0d0
	    sumby = 0.0d0
	    sumbz = 0.0d0
            sumfx = 0.0d0
	    sumfy = 0.0d0
	    sumfz = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|

            VoluReg = (0.5d0/3.0d0)*AreaCell(i)*(dsigv(k-1)+dsigv(k)) 

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              -------------------------------
!              Function at the points
               funv1T = funv(jv1,k)
               funv2T = funv(jv2,k)
!              _______________________________
!             |                               |
!             |      Region faces: 1,2,3      |
!             |_______________________________|

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
               sumax = sumax + (1.0d0/3.0d0)*nxAreaoVol
               sumay = sumay + (1.0d0/3.0d0)*nyAreaoVol
               sumaz = sumaz + (1.0d0/3.0d0)*nzAreaoVol
!              ------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               (funv1T+funv2T)
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               (funv1T+funv2T)
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               (funv1T+funv2T) 
!              _______________________________
!             |                               |
!             |      Region faces: 4,5,6      |
!             |_______________________________|

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k+1)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k+1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
               sumbx = sumbx + (1.0d0/3.0d0)*nxAreaoVol
               sumby = sumby + (1.0d0/3.0d0)*nyAreaoVol
               sumbz = sumbz + (1.0d0/3.0d0)*nzAreaoVol
!              ------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               (funv1T+funv2T)
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               (funv1T+funv2T)
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               (funv1T+funv2T) 
            enddo

!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|

!           ----------------------------
!           Constant of cell center i
            ax(i,k,4) = sumax
	    ay(i,k,4) = sumay
            az(i,k,4) = sumaz
!           -------------------------------
!           Constant of neighbor X^T
            bx(i,k,4) = sumbx
	    by(i,k,4) = sumby
            bz(i,k,4) = sumbz
!           ----------------------------
!           Known vertex values      
            fx(i,k,4) = sumfx
	    fy(i,k,4) = sumfy
            fz(i,k,4) = sumfz

!           ____________________________________________________
!          |                                                    |
!          |             Vertical BOTTOM neighbor               |
!          |____________________________________________________|

            sumax = 0.0d0
	    sumay = 0.0d0
	    sumaz = 0.0d0
            sumbx = 0.0d0
	    sumby = 0.0d0
	    sumbz = 0.0d0
            sumfx = 0.0d0
	    sumfy = 0.0d0
	    sumfz = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|

            VoluReg = (0.5d0/3.0d0)*AreaCell(i)*(dsigv(k-1)+dsigv(k)) 

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              -------------------------------
!              Function at the points
               funv1B = funv(jv1,k-1)
               funv2B = funv(jv2,k-1)
!              _______________________________
!             |                               |
!             |      Region faces: 1,2,3      |
!             |_______________________________|

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k-1)-sig(k-1)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k-1)-sig(k-1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
               sumbx = sumbx + (1.0d0/3.0d0)*nxAreaoVol
               sumby = sumby + (1.0d0/3.0d0)*nyAreaoVol
               sumbz = sumbz + (1.0d0/3.0d0)*nzAreaoVol
!              ------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               (funv1B+funv2B)
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               (funv1B+funv2B)
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               (funv1B+funv2B) 
!              _______________________________
!             |                               |
!             |      Region faces: 4,5,6      |
!             |_______________________________|

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k-1)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k-1)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
               sumax = sumax + (1.0d0/3.0d0)*nxAreaoVol
               sumay = sumay + (1.0d0/3.0d0)*nyAreaoVol
               sumaz = sumaz + (1.0d0/3.0d0)*nzAreaoVol
!              ------------------------------
!              Known vertex values              
	       sumfx = sumfx + (1.0d0/3.0d0)*nxAreaoVol*&
                               (funv1B+funv2B)
	       sumfy = sumfy + (1.0d0/3.0d0)*nyAreaoVol*&
                               (funv1B+funv2B)
	       sumfz = sumfz + (1.0d0/3.0d0)*nzAreaoVol*&
                               (funv1B+funv2B) 
            enddo

!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|

!           ----------------------------
!           Constant of cell center i
            ax(i,k,5) = sumax
	    ay(i,k,5) = sumay
            az(i,k,5) = sumaz
!           -------------------------------
!           Constant of neighbor X^T
            bx(i,k,5) = sumbx
	    by(i,k,5) = sumby
            bz(i,k,5) = sumbz
!           ----------------------------
!           Known vertex values      
            fx(i,k,5) = sumfx
	    fy(i,k,5) = sumfy
            fz(i,k,5) = sumfz
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientEdgeReg2'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
