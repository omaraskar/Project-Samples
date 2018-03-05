!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               CALCULATION OF THE GRADIENT AT EACH FACE              !
!                              March 2013                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientEdge3D(ax,ay,az,bx,by,bz,fx,fy,fz,  &
                                 fun,xc,yc,sig,dsig,No_cp,    &
                                 funv,xv,yv,sigv,dsigv,No_vp)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program approximates the gradient of a function fun at      !
!    the 5 faces of the prism (horizontal: j=1,2,3, top: j=4 and      !
!    bottom: j=5) taking an approximations of the integral over a     !
!    3D region surronding each face (Appendix B).                     !
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
!  | --> fun     |(N_CELL,NZ)| Function "fun" at the center element|  !
!  | --> funv    |(N_VERT,NZ)| Function "fun" at the vertices      |  !
!  | --> No_cp   |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  | --> No_vp   |(N_VERT,3) | Node No. of surrounding three cell  |  !
!  | --> sig     |(NZ)       | Vertical location of the cell center|  !
!  | --> sigv    |(NZ-1)     | Vertical location of the vertices   |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k+1)               |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    N_VERT  |  Number of the vertices in the horizontal domain |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |____________|__________________________________________________|  !
!  |  areaTriaH |(N_CELL0,3,NZglobal,8)| area of the triangle      |  !
!  |  nxTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle x  |  !
!  |  nyTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle y  |  !
!  |  nzTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle z  |  !
!  |  areaTriaT |(N_CELL0,NZglobal,6)  | area of the top triangles |  !
!  |  nxTriaT   |(N_CELL0,NZglobal,6)  | normal x top triangle     |  !
!  |  nyTriaT   |(N_CELL0,NZglobal,6)  | normal y top triangle     |  !
!  |  nzTriaT   |(N_CELL0,NZglobal,6)  | normal z top triangle     |  !
!  |  areaTriaB |(N_CELL0,NZglobal,6)  | area of the bottom trian. |  !
!  |  nxTriaB   |(N_CELL0,NZglobal,6)  | normal x bottom triangle  |  !
!  |  nyTriaB   |(N_CELL0,NZglobal,6)  | normal y bottom triangle  |  !
!  |  nzTriaB   |(N_CELL0,NZglobal,6)  | normal z bottom triangle  |  !
!  |____________|______________________|___________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !  
!  | funvp1      | function at the triangle vertex point 1         |  !
!  | funvp2      | function at the triangle vertex point 2         |  !
!  | funvp3      | function at the triangle vertex point 3         |  !
!  | sumfx       | real variable to add gradient component x       |  !
!  | sumfy       | real variable to add gradient component y       |  !
!  | sumfz       | real variable to add gradient component z       |  !
!  | jv1,jv2,jv3 | Number index of the vertices                    |  !
!  | s           | Loop counter of the number of triangles         |  !
!  |_____________|_________________________________________________|  !   
!  | f0,f1,f2    |  Function at the last 3 points from the boundary|  !
!  | h1,h2       |  length dsig of last 2 boundary points          |  !
!  | deno        |  denominator of the finite difference formula   |  !
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
      real*8 :: AreaReg,VoluReg
      real*8 :: nxArea,nyArea,nzArea
      real*8 :: dzup,dzdown,c0up,c0down,c1up,c1down

      integer:: jc,s,ss
      integer:: jj,jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~
!#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientEdge3D '
         print*,' '
!#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Calculation of the gradient                     !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1
         do i=1,N_CELL0

!           ____________________________________________________
!          |                                                    |
!          |             Horizontal neighbors: j=1:3            |
!          |____________________________________________________|

            do j=1,3
!              _______________________________
!             |                               |
!             |         Index & values        |
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
               VoluReg  = AreaReg*dsigv(k)  
!              _______________________________
!             |                               |
!             |      Region faces: 1 & 2      |
!             |_______________________________|

               sumax = 0.0d0
	       sumay = 0.0d0
	       sumbx = 0.0d0
	       sumby = 0.0d0
               sumfx = 0.0d0
	       sumfy = 0.0d0
	       sumfz = 0.0d0

               do s=1,2
                  if (s.eq.1) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv1,k)
                     funvT = funv(jv1,k+1)
!                    ----------------------------
!                    normal*Area
                     nxArea = (yc(jc)-yv(jv1))*dsigv(k)
                     nyArea =-(xc(jc)-xv(jv1))*dsigv(k)
                  elseif (s.eq.2) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv2,k)
                     funvT = funv(jv2,k+1)
!                    ----------------------------
!                    normal*Area
	             nxArea = (yv(jv2)-yc(jc))*dsigv(k)
	             nyArea =-(xv(jv2)-xc(jc))*dsigv(k)
                  endif
!                 -------------------------------
!                 Constant cell neighbor j 
                  sumbx = sumbx + 0.25d0*nxArea/VoluReg
                  sumby = sumby + 0.25d0*nyArea/VoluReg
!                 -------------------------------
!                 Known vertex values              
	          sumfx = sumfx + 0.25d0*nxArea/VoluReg*(funvB+funvT)
	          sumfy = sumfy + 0.25d0*nyArea/VoluReg*(funvB+funvT)
               enddo !s

!              _______________________________
!             |                               |
!             |      Region faces: 3 & 4      |
!             |_______________________________|

               do s=3,4
                  if (s.eq.3) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv2,k)
                     funvT = funv(jv2,k+1)
!                    ----------------------------
!                    normal*Area
	             nxArea =  (yc(i)-yv(jv2))*dsigv(k)
	             nyArea = -(xc(i)-xv(jv2))*dsigv(k)
                  elseif (s.eq.4) then 
!                    ----------------------------
!                    Function at the points
                     funvB = funv(jv1,k)
                     funvT = funv(jv1,k+1)
!                    ----------------------------
!                    normal*Area
	             nxArea =  (yv(jv1)-yc(i))*dsigv(k)
	             nyArea = -(xv(jv1)-xc(i))*dsigv(k)
                  endif
!                 -------------------------------
!                 Constant cell-center i 
                  sumax = sumax + 0.25d0*nxArea/VoluReg
                  sumay = sumay + 0.25d0*nyArea/VoluReg
!                 -------------------------------
!                 Known vertex values              
	          sumfx = sumfx + 0.25d0*nxArea/VoluReg*(funvB+funvT)
	          sumfy = sumfy + 0.25d0*nyArea/VoluReg*(funvB+funvT)
               enddo !s
!              _______________________________
!             |                               |
!             | Region face 5,6: Top + bottom |
!             |_______________________________|

!              ----------------------------
!              Function at the points
               funv1T = funv(jv1,k+1)
               funv2T = funv(jv2,k+1)
               funv1B = funv(jv1,k)
               funv2B = funv(jv2,k)
!             -------------------------------
!              Known vertex values              
               sumfz = sumfz + 0.25d0/dsigv(k)*((funv1T+funv2T)&
                                               -(funv1B+funv2B))
!              _______________________________
!             |                               |
!             |      Final contributions      |
!             |_______________________________|

!              ----------------------------
!              Constant of cell center i
               ax(i,k,j) = sumax
	       ay(i,k,j) = sumay
               az(i,k,j) = sumaz
!              -------------------------------
!              Constant of the neighbor j 
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

!           --------------------------------
!           Interpolation values top face
            dzup   = abs(sigv(k+2)-sig(k+1))
            dzdown = abs(sig(k+1)-sigv(k+1))
            c0up   = dzup/(dzdown+dzup)
            c0down = dzdown/(dzdown+dzup)

!           --------------------------------	
!           Interpolation values bottom face
            dzup   = abs(sigv(k+1)-sig(k))
            dzdown = abs(sig(k)-sigv(k))
            c1up   = dzup/(dzdown+dzup)
            c1down = dzdown/(dzdown+dzup)

!           _______________________________
!          |                               |
!          |         Faces : 1,2,3         |
!          |_______________________________|

            sumfx = 0.0d0
	    sumfy = 0.0d0
	    sumfz = 0.0d0

            do s=1,3
!              --------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss)
!              --------------------------------
!              Volume of the region 
               AreaReg = areaCell(i)             
               VoluReg = AreaReg*dsig(k)  
!              --------------------------------
!              normal*Area
               nxArea = (yv(jv2)-yv(jv1))*dsig(k)
               nyArea =-(xv(jv2)-xv(jv1))*dsig(k)
!              --------------------------------
!              Function at the top points
               funv1T = c0down*funv(jv1,k+2) + c0up*funv(jv1,k+1)
               funv2T = c0down*funv(jv2,k+2) + c0up*funv(jv2,k+1) 
!              --------------------------------
!              Function at the bottom points
               funv1B = c1down*funv(jv1,k+1) + c1up*funv(jv1,k)  
               funv2B = c1down*funv(jv2,k+1) + c1up*funv(jv2,k)  
!              --------------------------------
!              Known vertex values              
               sumfx = sumfx + 0.25d0*nxArea/VoluReg* &
                               (funv1T+funv2T+funv1B+funv2B)
	       sumfy = sumfy + 0.25d0*nyArea/VoluReg* &
                               (funv1T+funv2T+funv1B+funv2B)
            enddo !s

!           _______________________________
!          |                               |
!          | Region face 4,5: Top + Bottom |
!          |_______________________________|

!           ----------------------------
!           Function at the top points
            funv1T = c0down*funv(1,k+2) + c0up*funv(1,k+1)  
            funv2T = c0down*funv(2,k+2) + c0up*funv(2,k+1)
            funv3T = c0down*funv(3,k+2) + c0up*funv(3,k+1) 
!           ----------------------------
!           Function at the bottom points
            funv1B = c1down*funv(1,k+1) + c1up*funv(1,k) 
            funv2B = c1down*funv(2,k+1) + c1up*funv(2,k) 
            funv3B = c1down*funv(3,k+1) + c1up*funv(3,k) 
!           ----------------------------
!           Known vertex values              
            sumfz= sumfz+(1.0d0/3.0d0)/dsig(k)*((funv1T+funv2T+funv3T)&
                                               -(funv1B+funv2B+funv3B))
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|

!           ----------------------------
!           Known vertex values      
            fx(i,k,4) = sumfx
	    fy(i,k,4) = sumfy
            fz(i,k,4) = sumfz

!           ____________________________________________________
!          |                                                    |
!          |             Vertical BOTTOM neighbor               |
!          |____________________________________________________|

!           --------------------------------
!           Interpolation values top face
            dzup   = abs(sigv(k+1)-sig(k))
            dzdown = abs(sig(k)-sigv(k))
            c0up   = dzup/(dzdown+dzup)
            c0down = dzdown/(dzdown+dzup)

!           --------------------------------	
!           Interpolation values bottom face
            dzup   = abs(sigv(k)-sig(k-1))
            dzdown = abs(sig(k-1)-sigv(k-1))
            c1up   = dzup/(dzdown+dzup)
            c1down = dzdown/(dzdown+dzup)

!           _______________________________
!          |                               |
!          |         Faces : 1,2,3         |
!          |_______________________________|

            sumfx = 0.0d0
	    sumfy = 0.0d0
	    sumfz = 0.0d0

            do s=1,3
!              --------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss)
!              --------------------------------
!              Volume of the region 
               AreaReg = areaCell(i)             
               VoluReg = AreaReg*dsig(k-1)  
!              --------------------------------
!              normal*Area
               nxArea = (yv(jv2)-yv(jv1))*dsig(k)
               nyArea =-(xv(jv2)-xv(jv1))*dsig(k)
!              --------------------------------
!              Function at the top points
               funv1T = c0down*funv(jv1,k+1) + c0up*funv(jv1,k) 
               funv2T = c0down*funv(jv2,k+1) + c0up*funv(jv2,k)
!              --------------------------------
!              Function at the bottom points
               funv1B = c1down*funv(jv1,k) + c1up*funv(jv1,k-1) 
               funv2B = c1down*funv(jv2,k) + c1up*funv(jv2,k-1) 
!              --------------------------------
!              Known vertex values              
               sumfx = sumfx + 0.25d0*nxArea/VoluReg* &
                               (funv1B+funv1T+funv2B+funv2T)
	       sumfy = sumfy + 0.25d0*nyArea/VoluReg* &
                               (funv1B+funv1T+funv2B+funv2T)
            enddo !s

!           _______________________________
!          |                               |
!          | Region face 4,5: Top + Bottom |
!          |_______________________________|

!           ----------------------------
!           Function at the top points
            funv1T = c0down*funv(1,k+1) + c0up*funv(1,k) 
            funv2T = c0down*funv(2,k+1) + c0up*funv(2,k)
            funv3T = c0down*funv(3,k+1) + c0up*funv(3,k) 
!           ----------------------------
!           Function at the bottom points
            funv1B = c1down*funv(1,k) + c1up*funv(1,k-1) 
            funv2B = c1down*funv(2,k) + c1up*funv(2,k-1) 
            funv3B = c1down*funv(3,k) + c1up*funv(3,k-1) 
!           ----------------------------
!           normal*Area
            nzArea = AreaReg
!           ----------------------------
!           Known vertex values              
            sumfz= sumfz+(1.0d0/3.0d0)/dsig(k)*((funv1T+funv2T+funv3T)&
                                               -(funv1B+funv2B+funv3B))
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|

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

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientEdge3D '
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
