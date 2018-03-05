!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    TESTS: GRADIENT AT EDGE VALUES                   !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testGradientEdge(SolAppro,SolExact,SolError, &
                                  xc,yc,sig,dsig,No_cp,nbe,   &
                                  xv,yv,sigv,dsigv,No_vp,nbev)  
   
 
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine tests the approximation of the gradient at the   !
!    faces of the prism.                                              !
!                                                                     !
!---------------------------------------------------------------------!
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
!  | --> nbe     |(N_CELL,3) | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> funv    |(N_VERT,NZ)| Function "fun" at the vertices      |  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT,3) | Tag type of vertex points           |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Local varaibles:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !  
!  | funApprox   |(N_CELL0,NZ)| Approximation function             |  !
!  | funExact    |(N_CELL0,NZ)| Exact function                     |  !
!  | ErrorA1     |(N_CELL0,NZ)| Absolute error at each point       |  !
!  | ErrorR1     |(N_CELLO,NZ)| Relative error at each point       |  !
!  | MaxErrorA1  | (1:5,1:3)  | Maximum absolute error             |  !
!  | MaxErrorR1  | (1:5,1:3)  | Maximum realtive error             |  !
!  | DisplayThis | integer    | Tag to display the solution        |  !
!  |_____________|____________|____________________________________|  !
!                                                                     !
!     TecplotOption    = 0 Save the approximate solution              !
!                      = 1 Save the exact solution in tecplot files   !
!                      = 2 Save the errors in tecplot files           !
!                                                                     !
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

      real*8,dimension(:,:) :: SolAppro(N_CELL,NZ)
      real*8,dimension(:,:) :: SolExact(N_CELL,NZ)
      real*8,dimension(:,:) :: SolError(N_CELL,NZ)

      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 

      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 

!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________| 

      real*8,dimension(:),allocatable  :: phi2D
      real*8,dimension(:),allocatable  :: phiv2D
!     ----------------------------------------
      real*8,dimension(:,:),allocatable :: phi
      real*8,dimension(:,:),allocatable :: phiv
      real*8,dimension(:,:),allocatable :: funApprox
      real*8,dimension(:,:),allocatable :: funExact
      real*8,dimension(:,:),allocatable :: ErrorA1
      real*8,dimension(:,:),allocatable :: ErrorR1
      real*8,dimension(:,:) :: MaxErrorA1(1:5,1:3)
      real*8,dimension(:,:) :: MaxErrorR1(1:5,1:3)
      real*8,dimension(:,:) :: sumErrorA(1:5,1:3)
      real*8,dimension(:,:) :: sumErrorR(1:5,1:3)
!     ----------------------------------------
      real*8,dimension(:,:),allocatable :: funEdge1
      real*8,dimension(:,:),allocatable :: funEdge2
      real*8,dimension(:,:),allocatable :: funEdge3
!     ----------------------------------------
      real*8,dimension(:,:,:),allocatable :: dphidxE,dphidyE,dphidzE
      real*8,dimension(:,:,:),allocatable :: dphidxA,dphidyA,dphidzA
!     ----------------------------------------
      real*8 ,dimension(:,:,:),allocatable :: ax,ay,az
      real*8 ,dimension(:,:,:),allocatable :: bx,by,bz
      real*8 ,dimension(:,:,:),allocatable :: fx,fy,fz
!     ----------------------------------------
      real*8 :: dzT,dzB,dz0,ccT,ccB,cc0T,cc0B
      real*8 :: axB,axT,ax0,bxB,bxT,bx0
      real*8 :: ayB,ayT,ay0,byB,byT,by0
      real*8 :: azB,azT,az0,bzB,bzT,bz0
!     ----------------------------------------
      integer, parameter :: nx = NN
      integer, parameter :: ny = NN
      integer:: ii,jj
!     ----------------------------------------
      integer:: jv1,jv2,jv3
      integer:: jc,m,s,Nface,elem
!     ----------------------------------------
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: Neumanndfdn2D,Neumanndfdn3D
      real*8 :: x,y,z,fB,dfBdn,nnx,nny,nnz
!     ----------------------------------------
      real*8 :: dfdxExam2D,dfdyExam2D
      real*8 :: dfdxExam3D,dfdyExam3D,dfdzExam3D
!     ----------------------------------------
      integer,parameter :: RegionType    = 2
      integer,parameter :: UseExactVert  = 1
      integer,parameter :: TecplotOption = 2
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: Test gradientEdge'
      print*,' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(phi2D(N_CELL),phiv2D(N_VERT))

      allocate(phi(N_CELL,NZ), &
               phiv(N_VERT,NZ-1))

      allocate(funApprox(N_CELL0,NZ), &
               funExact(N_CELL0,NZ),  &
               ErrorA1(N_CELL0,NZ),   &
               ErrorR1(N_CELL0,NZ))

      allocate(dphidxE(N_CELL0,NZ,5), &
               dphidyE(N_CELL0,NZ,5), &
               dphidzE(N_CELL0,NZ,5), &
               dphidxA(N_CELL0,NZ,5), &
               dphidyA(N_CELL0,NZ,5), &
               dphidzA(N_CELL0,NZ,5))

      allocate(ax(N_CELL0,NZ,5),  &
               ay(N_CELL0,NZ,5),  &
               az(N_CELL0,NZ,5),  &
               bx(N_CELL0,NZ,5),  &
               by(N_CELL0,NZ,5),  &
               bz(N_CELL0,NZ,5),  &               
               fx(N_CELL0,NZ,5),  &
               fy(N_CELL0,NZ,5),  &
               fz(N_CELL0,NZ,5))  

      allocate(funEdge1(N_CELL0,NZ), &
               funEdge2(N_CELL0,NZ), &
               funEdge3(N_CELL0,NZ))

!*********************************************************************!
!                                                                     !
!                               Test 2D                               !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.2) THEN
!         _____________________________________________________
!        |                                                     |
!        |            2D Function at cell-centers              |
!        |_____________________________________________________|

!        ----------------------------------------------
!        Function: inside points
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
            phi2D(i) = funSolExam2D(x,y)
         enddo
!        ----------------------------------------------
!        Boundary Condition of the cell-center points    

         call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!        ----------------------------------------------
!        3D format     
         do k=1,NZ  
            do i=1,N_CELL
               phi(i,k) = phi2D(i)
            enddo
         enddo
!         _____________________________________________________
!        |                                                     |
!        |           Function phiv at the vertices             |  
!        |_____________________________________________________|

!        ______________________________________________________
!        Using the exact function
         if (UseExactVert.eq.1) then
            do nv=1,N_VERT
               x = xv(nv)
               y = yv(nv)
               phiv2D(nv) = funSolExam2D(x,y)
            enddo
!        ______________________________________________________
!        Using vertex interpolation
         else
!           ------------------------------------------------
!           Interpolation of the inside vertex points 
            call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                                 phi2D,xc,yc,No_cp,nbe)
!           ------------------------------------------------
!           Boundary Conditions of the vertex points 
            call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                            phi2D,xc,yc,No_cp,nbe)
         endif

!        ______________________________________________________
!        3D format     
         do k=1,NZ-1  
            do nv=1,N_VERT
               phiv(nv,k) = phiv2D(nv)
            enddo
         enddo
!         _____________________________________________________
!        |                                                     |
!        |           Exact gradient at the point (xe,ye)       |  
!        |_____________________________________________________|

         DO k=1,NZ  
            do i=1,N_CELL0
               do j=1,3
                  x = xme(i,j)
                  y = yme(i,j)
!                 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  if (nbe(i).gt.0) then
	              jc=No_cp(i,j)
	              if (jc.gt.N_CELL0) then
                         x = xe(i,j)
                         y = ye(i,j)
                      endif
                  endif
!                 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  dphidxE(i,k,j) = dfdxExam2D(x,y)
                  dphidyE(i,k,j) = dfdyExam2D(x,y)
	       enddo
            enddo
         ENDDO
!         _____________________________________________________
!        |                                                     |
!        |         Coefficients of the gradient of phi         |  
!        |_____________________________________________________|

         call grandientEdge2D(ax,ay,az,bx,by,bz,fx,fy,fz,     &
                              phi,xc,yc,sig,dsig,No_cp,       &
                              phiv,xv,yv,sigv,dsigv,No_vp)     

!         _____________________________________________________
!        |                                                     |
!        |          Approximation of the gradient              |  
!        |_____________________________________________________|

         DO k=2,NZ-1  
            do i=1,N_CELL0         
               do j=1,3
                  jc = No_cp(i,j)
                  dphidxA(i,k,j) = fx(i,k,j)&
                                  +ax(i,k,j)*phi(i,k) &
                                  +bx(i,k,j)*phi(jc,k) 
                  dphidyA(i,k,j) = fy(i,k,j)&
                                  +ay(i,k,j)*phi(i,k)&
                                  +by(i,k,j)*phi(jc,k)                   
	       enddo
            enddo
         ENDDO
!         _____________________________________________________
!        |                                                     |
!        |          BC of the gradient at the edge             |  
!        |_____________________________________________________|

          IF (ChooseBoundary.eq.2) THEN
             do i=1,N_CELL0
                if (nbe(i).ne.0) then	
   	           do j=1,3
	              nc = No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                         x = xe(i,j)
                         y = ye(i,j)
                         do k=1,NZ
                            dphidxA(i,k,j) = dfdxExam2D(x,y)
                            dphidyA(i,k,j) = dfdyExam2D(x,y)
                         enddo
                      endif 
	           enddo
                endif
             enddo
          ENDIF

!*********************************************************************!
!                                                                     !
!                               Test 3D                               !
!                                                                     !
!*********************************************************************!

      ELSEIF (TestDimension.eq.3) THEN
!         _____________________________________________________
!        |                                                     |
!        |             Function at cell-centers                |
!        |_____________________________________________________|

!        ----------------------------------------------
!        Function: inside points
         do i=1,N_CELL0
            do k=1,NZ 
               x = xc(i)
               y = yc(i)
               z = sig(k)
               phi(i,k) = funSolExam3D(x,y,z) 
            enddo
         enddo
!        ----------------------------------------------
!        Boundary Condition of the cell-center points  

         call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe)

!         _____________________________________________________
!        |                                                     |
!        |           Function phiv at the vertices             |  
!        |_____________________________________________________|

!        ______________________________________________________
!        Using the exact function
         if (UseExactVert.eq.1) then
            do nv=1,N_VERT
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)
                  phiv(nv,k) = funSolExam3D(x,y,z)
               enddo
            enddo
!        ______________________________________________________
!        Using vertex interpolation
         else
!           ------------------------------------------------
!           Interpolation of the inside vertex points 
            call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                                 phi,xc,yc,sig,dsig,No_cp,nbe)
!           ------------------------------------------------
!           Boundary Conditions of the vertex points 
            call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                            phi,xc,yc,sig,dsig,No_cp,nbe)
         endif
!         _____________________________________________________
!        |                                                     |
!        |           Exact gradient at the point (xe,ye)       |  
!        |_____________________________________________________|

         DO k=2,NZ-1  
            do i=1,N_CELL0
!              ---------------------------------------	
!              Horizontal             
               do j=1,3
                  x = xme(i,j)
                  y = yme(i,j)
!                 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  if (nbe(i).gt.0) then
	              jc=No_cp(i,j)
	              if (jc.gt.N_CELL0) then
                         x = xe(i,j)
                         y = ye(i,j)
                      endif
                  endif
!                 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                  z = sig(k)
                  dphidxE(i,k,j) = dfdxExam3D(x,y,z)
                  dphidyE(i,k,j) = dfdyExam3D(x,y,z) 
                  dphidzE(i,k,j) = dfdzExam3D(x,y,z) 
	       enddo
!              ---------------------------------------	
!              Vertical top
               x = xc(i)
               y = yc(i)
               z = sigv(k)
               dphidxE(i,k,4) = dfdxExam3D(x,y,z) 
               dphidyE(i,k,4) = dfdyExam3D(x,y,z) 
               dphidzE(i,k,4) = dfdzExam3D(x,y,z)
!              ---------------------------------------	
!              Vertical Bottom
               x = xc(i)
               y = yc(i)
               z = sigv(k-1)
               dphidxE(i,k,5) = dfdxExam3D(x,y,z) 
               dphidyE(i,k,5) = dfdyExam3D(x,y,z) 
               dphidzE(i,k,5) = dfdzExam3D(x,y,z) 
            enddo
         ENDDO

!         _____________________________________________________
!        |                                                     |
!        |          Approximation of the gradient Reg 1        |  
!        |_____________________________________________________|

         IF (RegionType.eq.1) THEN

            call grandientEdge3DReg1(ax,ay,az,bx,by,bz,fx,fy,fz, &
                                     phi,xc,yc,sig,dsig,No_cp,   &
                                     phiv,xv,yv,sigv,dsigv,No_vp)  
            DO k=2,NZ-1  
               do i=1,N_CELL0
!                 -----------------------------------
!                 Interpolation values
                  dzT = abs(sig(k+1)-sigv(k))
                  dz0 = abs(sigv(k)-sig(k))
                  dzB = abs(sigv(k-1)-sig(k-1))
                  ccT = dz0/(dz0+dzT)
                  cc0T= dzT/(dz0+dzT)
                  cc0B= dz0/(dz0+dzB)
                  ccB = dzB/(dz0+dzB)
!                  _______________________________
!                 |                               |   
!                 |     Horizontal interfaces     |
!                 |_______________________________|
             
                  do j=1,3
                     jc = No_cp(i,j)
!                    --------------------------------
!                    d(phi)/dx
                     axB = ax(i,k,j)*ccB        *phi(i,k-1)
                     ax0 = ax(i,k,j)*(cc0T+cc0B)*phi(i,k)
                     axT = ax(i,k,j)*ccT        *phi(i,k+1)
                     bxB = bx(i,k,j)*ccB        *phi(jc,k-1)
                     bx0 = bx(i,k,j)*(cc0T+cc0B)*phi(jc,k)
                     bxT = bx(i,k,j)*ccT        *phi(jc,k+1)
                     dphidxA(i,k,j) = fx(i,k,j)   &   
                                     +axB+ax0+axT &
                                     +bxB+bx0+bxT                                    
!                    --------------------------------
!                    d(phi)/dy
                     ayB = ay(i,k,j)*ccB        *phi(i,k-1)
                     ay0 = ay(i,k,j)*(cc0T+cc0B)*phi(i,k)
                     ayT = ay(i,k,j)*ccT        *phi(i,k+1)
                     byB = by(i,k,j)*ccB        *phi(jc,k-1)
                     by0 = by(i,k,j)*(cc0T+cc0B)*phi(jc,k)
                     byT = by(i,k,j)*ccT        *phi(jc,k+1)
                     dphidyA(i,k,j) = fy(i,k,j)   &
                                     +ayB+ay0+ayT &
                                     +byB+by0+byT 
!                    --------------------------------
!                    d(phi)/dz
                     azB = 0.25d0/dsigv(k)*(-ccB)     *phi(i,k-1)
                     az0 = 0.25d0/dsigv(k)*(cc0T-cc0B)*phi(i,k)
                     azT = 0.25d0/dsigv(k)*ccT        *phi(i,k+1)
                     bzB = 0.25d0/dsigv(k)*(-ccB)     *phi(jc,k-1)
                     bz0 = 0.25d0/dsigv(k)*(cc0T-cc0B)*phi(jc,k)
                     bzT = 0.25d0/dsigv(k)*ccT        *phi(jc,k+1)
   
                     dphidzA(i,k,j) = fz(i,k,j)   &
                                     +azB+az0+azT &
                                     +bzB+bz0+bzT 
	          enddo
!                  _______________________________
!                 |                               |
!                 |     Vertical top interface    |
!                 |_______________________________|

                  dphidxA(i,k,4) = fx(i,k,4)
                  dphidyA(i,k,4) = fy(i,k,4)                             
                  dphidzA(i,k,4) = fz(i,k,4)
!                  _______________________________
!                 |                               |
!                 |   Vertical bottom interface   |
!                 |_______________________________|

                  dphidxA(i,k,5) = fx(i,k,5)
                  dphidyA(i,k,5) = fy(i,k,5) 
                  dphidzA(i,k,5) = fz(i,k,5) 
               enddo
            ENDDO
!         _____________________________________________________
!        |                                                     |
!        |          Approximation of the gradient Reg 2        |  
!        |_____________________________________________________|

         ELSEIF (RegionType.eq.2) THEN

            call grandientEdge3DReg2(ax,ay,az,bx,by,bz,fx,fy,fz,  &
                                     phi,xc,yc,sig,dsig,No_cp,    &
                                     phiv,xv,yv,sigv,dsigv,No_vp)  

            DO k=2,NZ-1  
               do i=1,N_CELL0
!                  _______________________________
!                 |                               |
!                 |     Horizontal interfaces     |
!                 |_______________________________|
             
                  do j=1,3
                     jc = No_cp(i,j)
!                    ------------------------------
!                    d(phi)/dx
                     dphidxA(i,k,j) = fx(i,k,j)          &
                                     +ax(i,k,j)*phi(i,k) &
                                     +bx(i,k,j)*phi(jc,k)                                   
!                    ------------------------------
!                    d(phi)/dy
                     dphidyA(i,k,j) = fy(i,k,j)          &
                                     +ay(i,k,j)*phi(i,k) &
                                     +by(i,k,j)*phi(jc,k) 
!                    ------------------------------
!                    d(phi)/dz
                     dphidzA(i,k,j) = fz(i,k,j)          &
                                     +az(i,k,j)*phi(i,k) &
                                     +bz(i,k,j)*phi(jc,k) 
	          enddo
!                  _______________________________
!                 |                               |
!                 |     Vertical top interface    |
!                 |_______________________________|

!                 ------------------------------
!                 d(phi)/dx
                  dphidxA(i,k,4) = fx(i,k,4)            &
                                  +ax(i,k,4)*phi(i,k)   &
                                  +bx(i,k,4)*phi(i,k+1)  
!                 ------------------------------
!                 d(phi)/dy
                  dphidyA(i,k,4) = fy(i,k,4)            &
                                  +ay(i,k,4)*phi(i,k)   &
                                  +by(i,k,4)*phi(i,k+1)  
!                 ------------------------------
!                 d(phi)/dz                            
                  dphidzA(i,k,4) = fz(i,k,4)            &
                                  +az(i,k,4)*phi(i,k)   &
                                  +bz(i,k,4)*phi(i,k+1) 
!                  _______________________________
!                 |                               |
!                 |   Vertical bottom interface   |
!                 |_______________________________|

!                 ------------------------------
!                 d(phi)/dx
                  dphidxA(i,k,5) = fx(i,k,5)            &
                                  +ax(i,k,5)*phi(i,k)   &
                                  +bx(i,k,5)*phi(i,k-1)  
!                 ------------------------------
!                 d(phi)/dy
                  dphidyA(i,k,5) = fy(i,k,5)            &
                                  +ay(i,k,5)*phi(i,k)   &
                                  +by(i,k,5)*phi(i,k-1)
!                 ------------------------------
!                 d(phi)/dz                              
                  dphidzA(i,k,5) = fz(i,k,5)            &
                                  +az(i,k,5)*phi(i,k)   &
                                  +bz(i,k,5)*phi(i,k-1) 
               enddo
            ENDDO
         ENDIF

      ENDIF

!*********************************************************************!
!                                                                     !
!                      Vertical boundary values                       !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |  We only need to update the central values from k=2    |
!     |  to k=NZ-1, the other values are calculated using the  |
!     |  corresponding boundary conditions. In this case we    |
!     |  I will assign the correct values at those points.     |
!     |________________________________________________________|

      k=1  
      do i=1,N_CELL0
         do j=1,5
            dphidxE(i,k,j) = 0.0d0
            dphidyE(i,k,j) = 0.0d0 
            dphidzE(i,k,j) = 0.0d0 
            dphidxA(i,k,j) = dphidxE(i,k,j)
            dphidyA(i,k,j) = dphidyE(i,k,j) 
            dphidzA(i,k,j) = dphidzE(i,k,j) 
        enddo
      enddo
      k=NZ  
      do i=1,N_CELL0
         do j=1,5
            dphidxE(i,k,j) = 0.0d0
            dphidyE(i,k,j) = 0.0d0 
            dphidzE(i,k,j) = 0.0d0
            dphidxA(i,k,j) = dphidxE(i,k,j)
            dphidyA(i,k,j) = dphidyE(i,k,j) 
            dphidzA(i,k,j) = dphidzE(i,k,j) 
        enddo
      enddo
 
!*********************************************************************!
!                                                                     !
!                                 Error                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                Calculation of the errors               |
!     |________________________________________________________|

      if (TestDimension.eq.2) Nface = 3    
      if (TestDimension.eq.3) Nface = 5
      DO m=1,3
         Do s=1,Nface
!           --------------------------------------------------
!           Testing functions
            if (m.eq.1) then
               do k=1,NZ 
                  do i=1,N_CELL0 
                     funApprox(i,k) = dphidxA(i,k,s) 
                     funExact(i,k)  = dphidxE(i,k,s)
                  enddo
               enddo
            elseif (m.eq.2) then
               do k=1,NZ 
                  do i=1,N_CELL0 
                     funApprox(i,k) = dphidyA(i,k,s) 
                     funExact(i,k)  = dphidyE(i,k,s)
                  enddo
               enddo
            elseif (m.eq.3) then
               do k=1,NZ 
                  do i=1,N_CELL0 
                     funApprox(i,k) = dphidzA(i,k,s) 
                     funExact(i,k)  = dphidzE(i,k,s)
                  enddo
               enddo
            endif

!           --------------------------------------------------
!           Error (Maximum norm: p=infty) 
            maxErrorA1(s,m) = 0.0d0
            maxErrorR1(s,m) = 0.0d0
            sumErrorA(s,m) = 0.0d0
            sumErrorR(s,m) = 0.0d0
            do k=1,NZ 
               do i=1,N_CELL0 
!                 -----------------
!                 Absolute
                  ErrorA1(i,k)= abs(funApprox(i,k)-funExact(i,k))
                  maxErrorA1(s,m) = max(maxErrorA1(s,m),ErrorA1(i,k))
                  sumErrorA(s,m)  = sumErrorA(s,m) +(ErrorA1(i,k))**2
!                 -----------------
!                 Relative
                  ErrorR1(i,k) = abs(funExact(i,k)) 
                  maxErrorR1(s,m) = max(maxErrorR1(s,m),ErrorR1(i,k))
                  sumErrorR(s,m)  = sumErrorR(s,m) +(ErrorR1(i,k))**2
               enddo
            enddo
            maxErrorR1(s,m) = maxErrorA1(s,m)/maxErrorR1(s,m)
            sumErrorA(s,m)  = dsqrt(sumErrorA(s,m)/(N_CELL0*(NZ-1)))
            sumErrorR(s,m)  = dsqrt(sumErrorR(s,m)/(N_CELL0*(NZ-1)))
            sumErrorR(s,m)  = sumErrorA(s,m)/sumErrorR(s,m)
         Enddo
      ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                      Display solution                  |
!     |________________________________________________________|

9     format(t10,i3,t15,e10.3,t27,e10.3,t39,e10.3,t51,e10.3)
8     format(t10,60a)
      write(*,8)'===================================================='
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'                TEST: Gradient at edge              '
      write(*,'(t30,a4,i3)') ' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'Face  Max norm    L-2 norm    Max norm    L-2 norm  '
      m=1
      write(*,8)'----------------------------------------------------'
      write(*,8)'                        dphi/dx                     '
      write(*,8)'                                                    '
      do s=1,Nface
         write(*,9) s,maxErrorA1(s,m),sumErrorA(s,m),&
                      maxErrorR1(s,m),sumErrorR(s,m)
      enddo
      m=2
      write(*,8)'----------------------------------------------------'
      write(*,8)'                        dphi/dy                     '
      write(*,8)'                                                    '
      do s=1,Nface
         write(*,9) s,maxErrorA1(s,m),sumErrorA(s,m),&
                      maxErrorR1(s,m),sumErrorR(s,m)
      enddo
      if (TestDimension.eq.3) then
      m=3
      write(*,8)'----------------------------------------------------'
      write(*,8)'                        dphi/dz                     '
      write(*,8)'                                                    '
      do s=1,Nface
         write(*,9) s,maxErrorA1(s,m),sumErrorA(s,m),&
                      maxErrorR1(s,m),sumErrorR(s,m)
      enddo
      endif
      write(*,8)'----------------------------------------------------'
      write(*,8)'===================================================='
      write(*,*) ''

!      ________________________________________________________
!     |                                                        |
!     |       Solution to display in tecplot (cell-center)     |
!     |________________________________________________________|
 
      do k=1,NZ 
         do i=1,N_CELL0 
            SolAppro(i,k) = dphidxE(i,k,1)
            SolExact(i,k) = dphidxA(i,k,1)
            SolError(i,k) = abs(dphidxA(i,k,1)-dphidxE(i,k,1))
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |       Solution to display in tecplot (edges)           |
!     |________________________________________________________|

!     ---------------------------------
!     Exact solution in tecplot files        
      if (TecplotOption.eq.1) then
         do k=1,NZ 
            do i=1,N_CELL0 
               funEdge1(i,k) = dphidxE(i,k,1)
               funEdge2(i,k) = dphidxE(i,k,2)
               funEdge3(i,k) = dphidxE(i,k,3)
            enddo
         enddo
!     ---------------------------------
!     Error in tecplot files        
      elseif (TecplotOption.eq.2) then
         do k=1,NZ 
            do i=1,N_CELL0 
               funEdge1(i,k) = abs(dphidxA(i,k,1)-dphidxE(i,k,1))
               funEdge2(i,k) = abs(dphidxA(i,k,2)-dphidxE(i,k,2))
               funEdge3(i,k) = abs(dphidxA(i,k,3)-dphidxE(i,k,3))
            enddo
         enddo
!     ---------------------------------
!     Approximation in tecplot files 
      else
         do k=1,NZ 
            do i=1,N_CELL0 
               funEdge1(i,k) = dphidxA(i,k,1)
               funEdge2(i,k) = dphidxA(i,k,2)
               funEdge3(i,k) = dphidxA(i,k,3)
            enddo
         enddo
      endif


      call SavetecEdge(funEdge1,funEdge2,funEdge3,sig)


!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(phi2D,phiv2D)
      deallocate(phi,phiv)
      deallocate(ax,ay,az,bx,by,bz,fx,fy,fz)
      deallocate(dphidxE,dphidyE,dphidzE,& 
                 dphidxA,dphidyA,dphidzA)
      deallocate(funApprox,funExact,ErrorA1,ErrorR1)
      deallocate(funEdge1,funEdge2,funEdge3)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'),'<<<<< End   subroutine: Test gradientEdge'
      print*,' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
