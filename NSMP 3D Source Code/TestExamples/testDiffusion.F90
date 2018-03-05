!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           TEST DIFFUSION                            !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testDiffusion(SolAppro,SolExact,SolError,      &
                               SolApprov,SolExactv,SolErrorv,   & 
                               xc,yc,sig,dsig,No_cp,nbe,        &
                               xv,yv,sigv,dsigv,No_vp,nbev)              
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program test the diffusion term of a given example phi.     !  
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name      |    Size     | Description                     |  !  
!  |_______________|_____________|_________________________________|  !
!  | <-- SolAppro  |(N_CELL,NZ)  | Approximate solution cell-center|  !
!  | <-- SolExact  |(N_CELL,NZ)  | Exact solution cell-center      |  !  
!  | <-- SolError  |(N_CELL,NZ)  | Error cell-center               |  !    
!  | <-- SolApprov |(N_VERT,NZ-1)| Approximate solution vertex     |  !
!  | <-- SolExactv |(N_VERT,NZ-1)| Exact solution vertex           |  ! 
!  | <-- SolErrorv |(N_VERT,NZ-1)| Error vertex                    |  ! 
!  |_______________|_____________|_________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbev   |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!                                                                     !
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

      real*8,dimension(:,:) :: SolAppro(N_CELL,NZ)
      real*8,dimension(:,:) :: SolExact(N_CELL,NZ)
      real*8,dimension(:,:) :: SolError(N_CELL,NZ)
!     --------------------------------------
      real*8,dimension(:,:) :: SolApprov(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolExactv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolErrorv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     --------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:) :: phi2D(N_CELL)
      real*8,dimension(:) :: phiv2D(N_VERT)
      real*8,dimension(:) :: Am02D(N_CELL0)
      real*8,dimension(:) :: Am12D(N_CELL0)
      real*8,dimension(:) :: Am22D(N_CELL0)
      real*8,dimension(:) :: Am32D(N_CELL0)
      real*8,dimension(:) :: Bmv12D(N_CELL0) 
      real*8,dimension(:) :: Bmv22D(N_CELL0) 
      real*8,dimension(:) :: Bmv32D(N_CELL0) 
      real*8,dimension(:) :: Nm2D(N_CELL0)
      real*8,dimension(:) :: Gamx2D(N_CELL)
      real*8,dimension(:) :: Gamy2D(N_CELL)
      real*8,dimension(:) :: Gam2D(N_CELL)
!     ----------------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv(N_CELL0,NZ) 
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Nm(N_CELL0,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)  
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)  
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)  
      real*8,dimension(:,:) :: Gam(N_CELL,NZ)  
!     ----------------------------------------
      real*8,dimension(:,:) :: diffA(N_CELL0,NZ)
!     ----------------------------------------
      real*8 :: x,y,z
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: funDiffExam2D,funDiffExam3D
      real*8 :: dphidx,dphidy,dphidz,suma
      real*8 :: dfdxExam2D,dfdyExam2D
      real*8 :: dfdxExam3D,dfdyExam3D,dfdzExam3D
      real*8 :: Volume
!     ----------------------------------------
      real*8 ,dimension(:,:,:),allocatable :: ax,ay,az
      real*8 ,dimension(:,:,:),allocatable :: bx,by,bz
      real*8 ,dimension(:,:,:),allocatable :: fx,fy,fz
!     ----------------------------------------
      real*8,dimension(:,:) :: funApprox(N_CELL,NZ)
      real*8,dimension(:,:) :: funExact(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
!     ----------------------------------------
      integer:: jv1,jv2,jv3,jj
      integer:: jc1,jc2,jc3,jc
!     ----------------------------------------
      integer,parameter :: RegionType       = 2
      integer,parameter :: UseDiffMethod    = 1
      integer,parameter :: UseExactGradient = 0 
      integer,parameter :: UseExactVert     = 0

!     ---------------------------------------- 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,*) ''
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: TEST DIFFUSION'
      write(*,*) ''
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(ax(N_CELL0,NZ,5),  &
               ay(N_CELL0,NZ,5),  &
               az(N_CELL0,NZ,5),  &
               bx(N_CELL0,NZ,5),  &
               by(N_CELL0,NZ,5),  &
               bz(N_CELL0,NZ,5),  &               
               fx(N_CELL0,NZ,5),  &
               fy(N_CELL0,NZ,5),  &
               fz(N_CELL0,NZ,5))  

!*********************************************************************!
!                                                                     !
!                        ---   2D CASE  ---                           !
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
!        ________________________________________________________
!       |                                                        |
!       |   Approximation using directly the GradientEdge.F90    |
!       |________________________________________________________|


         if (UseDiffMethod.eq.0) then

!           _____________________________________________________
!           Gradient at the edge      

            call grandientEdge2D(ax,ay,az,bx,by,bz,fx,fy,fz,     &
                                 phi,xc,yc,sig,dsig,No_cp,       &
                                 phiv,xv,yv,sigv,dsigv,No_vp)    

!           _____________________________________________________
!           Solution             

            do k=2,NZ-1  
               do i=1,N_CELL0 
                  suma = 0.0d0
                  do j=1,3
!                    --------------------------------
!                    Cell-center index
                     jc = No_cp(i,j)
!                    --------------------------------
!                    Vertex index               
	             jj=j+1
	             if (jj.gt.3) jj=jj-3
	             jv1 = No_vp(i,j)
	             jv2 = No_vp(i,jj)
!                    --------------------------------
!                    Gradient      
                     if (UseExactGradient.eq.1) then
                        x = 0.5d0*(xv(jv1)+xv(jv2))
                        y = 0.5d0*(yv(jv1)+yv(jv2))
                        dphidx = dfdxExam2D(x,y)
                        dphidy = dfdyExam2D(x,y)
                     else
                        dphidx = fx(i,k,j)&
                                +ax(i,k,j)*phi2D(i) &
                                +bx(i,k,j)*phi2D(jc) 
                        dphidy = fy(i,k,j)&
                                +ay(i,k,j)*phi2D(i)&
                                +by(i,k,j)*phi2D(jc) 
                     endif
!                    --------------------------------
!                    Contribution 
                     suma = suma + dphidx*(yv(jv2)-yv(jv1)) &
                                  -dphidy*(xv(jv2)-xv(jv1))                  
	          enddo
!                 -----------------------------------
!                 Approximated solution
                  diffA(i,k) = suma
               enddo
            enddo

!        ________________________________________________________
!       |                                                        |
!       |        Approximation using  diffusion2D.F90            |
!       |________________________________________________________|

         elseif (UseDiffMethod.eq.1) then
!           _____________________________________________________
!           Initialization         
            do i=1,N_CELL0 	
               Am02D(i)  = 0.0d0 
               Am12D(i)  = 0.0d0
               Am22D(i)  = 0.0d0
               Am32D(i)  = 0.0d0
               Bmv12D(i) = 0.0d0 
               Bmv22D(i) = 0.0d0 
               Bmv32D(i) = 0.0d0 
            enddo
!           _____________________________________________________
!           Diffusive coefficients
            do i=1,N_CELL	
               Gamx2D(i) = 1.0d0 
               Gamy2D(i) = 1.0d0
            enddo
!           _____________________________________________________
!           Diffusion constants    
            call diffusion2D(Am02D,Am12D,Am22D,Am32D, &
                             Bmv12D,Bmv22D,Bmv32D,    &
                             Gamx2D,Gamy2D,           &
                             xc,yc,No_cp,nbe,         &
                             xv,yv,No_vp,nbev)
!           _____________________________________________________
!           Solution 
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              ----------------          
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               do k=1,NZ
                  diffA(i,k)=  Am02D(i)*phi2D(i)     &  
                             + Am12D(i)*phi2D(jc1)   &
                             + Am22D(i)*phi2D(jc2)   &
                             + Am32D(i)*phi2D(jc3)   &
!                            -----------------------
                             + Bmv12D(i)*phiv2D(jv1) &
                             + Bmv22D(i)*phiv2D(jv2) &
                             + Bmv32D(i)*phiv2D(jv3) 
               enddo
            enddo

!        ________________________________________________________
!       |                                                        |
!       |      Approximation using  diffusionNeumann2D.F90       |
!       |________________________________________________________|

         elseif (UseDiffMethod.eq.2) then
!           _____________________________________________________
!           Initialization         
            do i=1,N_CELL0 	
               Am02D(i)  = 0.0d0 
               Am12D(i)  = 0.0d0
               Am22D(i)  = 0.0d0
               Am32D(i)  = 0.0d0
               Bmv12D(i) = 0.0d0 
               Bmv22D(i) = 0.0d0 
               Bmv32D(i) = 0.0d0 
               Nm2D(i)   = 0.0d0 
            enddo
!           _____________________________________________________
!           Diffusive coefficients
            do i=1,N_CELL	
               Gamx2D(i) = 1.0d0 
               Gamy2D(i) = 1.0d0
            enddo
!           _____________________________________________________
!           Diffusion constants    
            if (ChooseBoundary.eq.2) then
               call diffusionNeumann2D(Am02D,Am12D,Am22D,Am32D,  &
                                       Bmv12D,Bmv22D,Bmv32D,Nm2D,&
                                       Gamx2D,Gamy2D,            &
                                       xc,yc,No_cp,nbe,          &
                                       xv,yv,No_vp,nbev)
            else
               print*,'Error!!!: This test is only for Neumann BC'
               stop
            endif
!           _____________________________________________________
!           Solution 
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              ----------------          
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               do k=1,NZ
                  diffA(i,k)=  Am02D(i)*phi2D(i)     &  
                             + Am12D(i)*phi2D(jc1)   &
                             + Am22D(i)*phi2D(jc2)   &
                             + Am32D(i)*phi2D(jc3)   &
!                            -----------------------
                             + Bmv12D(i)*phiv2D(jv1) &
                             + Bmv22D(i)*phiv2D(jv2) &
                             + Bmv32D(i)*phiv2D(jv3) &
!                            -----------------------
                             + Nm2D(i)
               enddo
            enddo

         endif !(UseDiffMethod)

!*********************************************************************!
!                                                                     !
!                        ---   3D CASE  ---                           !
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

!        ________________________________________________________
!       |                                                        |
!       |   Approximation using directly gradientEdgeReg2.F90    |
!       |________________________________________________________|

         if (UseDiffMethod.eq.0) then

!        _______________________________
!       |                               |
!       |    Region type 2: Piramids    |
!       |_______________________________|

         if (RegionType.eq.2) then

!           _____________________________________________________
!           Gradient at the edge  

            call grandientEdge3DReg2(ax,ay,az,bx,by,bz,fx,fy,fz, &
                                     phi,xc,yc,sig,dsig,No_cp,   &
                                     phiv,xv,yv,sigv,dsigv,No_vp)  
!           _____________________________________________________
!           Solution             

            DO k=2,NZ-1  
               do i=1,N_CELL0
                  suma = 0.0d0
!                 _______________________________
!                 Horizontal interfaces     
                  do j=1,3
!                    ----------------------------
                     jc = No_cp(i,j)
!                    ----------------------------
	             jj=j+1
	             if (jj.gt.3) jj=jj-3
	             jv1 = No_vp(i,j)
	             jv2 = No_vp(i,jj)
!                    ----------------------------
                     if (UseExactGradient.eq.1) then
                        x = 0.5d0*(xv(jv1)+xv(jv2))
                        y = 0.5d0*(yv(jv1)+yv(jv2))
                        z = sig(k)
                        dphidx = dfdxExam3D(x,y,z)
                        dphidy = dfdyExam3D(x,y,z) 
                        dphidz = dfdzExam3D(x,y,z) 
                     else
                        dphidx = ax(i,k,j)*phi(i,k)  &
                                +bx(i,k,j)*phi(jc,k) &
                                +fx(i,k,j)          
                        dphidy = ay(i,k,j)*phi(i,k)  &
                                +by(i,k,j)*phi(jc,k) &
                                +fy(i,k,j)
                        dphidz = az(i,k,j)*phi(i,k)  &
                                +bz(i,k,j)*phi(jc,k) &
                                +fz(i,k,j)
                     endif
!                    ----------------------------
                     suma = suma + dphidx*(yv(jv2)-yv(jv1))*dsigv(k) &
                                  -dphidy*(xv(jv2)-xv(jv1))*dsigv(k)  
	          enddo
!                 _______________________________
!                 Vertical top interface   
                  if (UseExactGradient.eq.1) then 
                     x = xc(i)
                     y = yc(i)
                     z = sigv(k)
                     dphidx = dfdxExam3D(x,y,z) 
                     dphidy = dfdyExam3D(x,y,z) 
                     dphidz = dfdzExam3D(x,y,z)
                  else                 
                     dphidx = ax(i,k,4)*phi(i,k)   &
                             +bx(i,k,4)*phi(i,k+1) & 
                             +fx(i,k,4)            
                     dphidy = ay(i,k,4)*phi(i,k)   &
                             +by(i,k,4)*phi(i,k+1) &  
                             +fy(i,k,4)            
                     dphidz = az(i,k,4)*phi(i,k)   &
                             +bz(i,k,4)*phi(i,k+1) &
                             +fz(i,k,4)
                  endif
                  suma = suma + dphidz*AreaCell(i)   
!                 _______________________________
!                 Vertical bottom interface   
                  if (UseExactGradient.eq.1) then 
                     x = xc(i)
                     y = yc(i)
                     z = sigv(k-1)
                     dphidx = dfdxExam3D(x,y,z) 
                     dphidy = dfdyExam3D(x,y,z) 
                     dphidz = dfdzExam3D(x,y,z)
                  else
                     dphidx = ax(i,k,5)*phi(i,k)   &
                             +bx(i,k,5)*phi(i,k-1) &
                             +fx(i,k,5)
                     dphidy = ay(i,k,5)*phi(i,k)   &
                             +by(i,k,5)*phi(i,k-1) &
                             +fy(i,k,5)
                     dphidz = az(i,k,5)*phi(i,k)   &
                             +bz(i,k,5)*phi(i,k-1) &
                             +fz(i,k,5)
                  endif
                  suma = suma - dphidz*AreaCell(i)
!                 _______________________________
!                 Approximated solution
                  diffA(i,k) = suma
               enddo
            ENDDO
            
         endif !(RegionType)

!        ________________________________________________________
!       |                                                        |
!       |           Approximation using diffusion3D.F90          |
!       |________________________________________________________|

         elseif (UseDiffMethod.eq.1) then

!           _____________________________________________________
!           Initialization 
       
            do k=1,NZ
               do i=1,N_CELL0 	
                  Am0(i,k)  = 0.0d0 
                  Am1(i,k)  = 0.0d0
                  Am2(i,k)  = 0.0d0
                  Am3(i,k)  = 0.0d0
                  AmT(i,k)  = 0.0d0
                  AmB(i,k)  = 0.0d0
                  Bmv1T(i,k)= 0.0d0 
                  Bmv2T(i,k)= 0.0d0 
                  Bmv3T(i,k)= 0.0d0 
                  Bmv1B(i,k)= 0.0d0 
                  Bmv2B(i,k)= 0.0d0 
                  Bmv3B(i,k)= 0.0d0 
               enddo
            enddo
!           ________________________________________________________
!           Diffusive coefficients
            do k=1,NZ
               do i=1,N_CELL	
                  Gamx(i,k) = 1.0d0 
                  Gamy(i,k) = 1.0d0 
                  Gamz(i,k) = 1.0d0 
               enddo
            enddo
!           _____________________________________________________
!           Diffusion constants

            call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                             Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                             Gamx,Gamy,Gamz,                      &
                             xc,yc,sig,dsig,No_cp,nbe,            &
                             xv,yv,sigv,dsigv,No_vp,nbev)   

!           _____________________________________________________
!           Solution (inside)

            do k=2,NZ-1
               do i=1,N_CELL0 	
                  jc1 = No_cp(i,1)
                  jc2 = No_cp(i,2)
                  jc3 = No_cp(i,3)
                  jv1 = No_vp(i,1)
                  jv2 = No_vp(i,2)
                  jv3 = No_vp(i,3)
                  diffA(i,k) =  (Am0(i,k)*phi(i,k)        &  
                               + Am1(i,k)*phi(jc1,k)      &
                               + Am2(i,k)*phi(jc2,k)      &
                               + Am3(i,k)*phi(jc3,k)      &
                               + AmT(i,k)*phi(i,k+1)      &       
                               + AmB(i,k)*phi(i,k-1))     &
!                              ----------------------------
                               + Bmv1T(i,k)*phiv(jv1,k)   &
                               + Bmv2T(i,k)*phiv(jv2,k)   &
                               + Bmv3T(i,k)*phiv(jv3,k)   &
                               + Bmv1B(i,k)*phiv(jv1,k-1) &
                               + Bmv2B(i,k)*phiv(jv2,k-1) &
                               + Bmv3B(i,k)*phiv(jv3,k-1)  
               enddo
            enddo
!        ________________________________________________________
!       |                                                        |
!       |      Approximation using  diffusionNeumann3D.F90       |
!       |________________________________________________________|

         elseif (UseDiffMethod.eq.2) then

!           _____________________________________________________
!           Initialization 
       
            do k=1,NZ
               do i=1,N_CELL0 	
                  Am0(i,k)  = 0.0d0 
                  Am1(i,k)  = 0.0d0
                  Am2(i,k)  = 0.0d0
                  Am3(i,k)  = 0.0d0
                  AmT(i,k)  = 0.0d0
                  AmB(i,k)  = 0.0d0
                  Bmv1T(i,k)= 0.0d0 
                  Bmv2T(i,k)= 0.0d0 
                  Bmv3T(i,k)= 0.0d0 
                  Bmv1B(i,k)= 0.0d0 
                  Bmv2B(i,k)= 0.0d0 
                  Bmv3B(i,k)= 0.0d0 
                  Nm(i,k)   = 0.0d0 
               enddo
            enddo
!           ________________________________________________________
!           Diffusive coefficients
            do k=1,NZ
               do i=1,N_CELL	
                  Gamx(i,k) = 1.0d0 
                  Gamy(i,k) = 1.0d0 
                  Gamz(i,k) = 1.0d0 
               enddo
            enddo
!           _____________________________________________________
!           Diffusion constants
            if (ChooseBoundary.eq.2) then
               call diffusionNeumann3D(Am0,Am1,Am2,Am3,AmT,AmB,    &
                                       Bmv1T,Bmv2T,Bmv3T,          &
                                       Bmv1B,Bmv2B,Bmv3B,Nm,       & 
                                       Gamx,Gamy,Gamz,             &
                                       xc,yc,sig,dsig,No_cp,nbe,   &
                                       xv,yv,sigv,dsigv,No_vp,nbev)   
            else
               print*,'Error!!!: This test is only for Neumann BC'
               stop
            endif
!           _____________________________________________________
!           Solution (inside)

            do k=2,NZ-1
               do i=1,N_CELL0 	
                  jc1 = No_cp(i,1)
                  jc2 = No_cp(i,2)
                  jc3 = No_cp(i,3)
                  jv1 = No_vp(i,1)
                  jv2 = No_vp(i,2)
                  jv3 = No_vp(i,3)
                  diffA(i,k) =  (Am0(i,k)*phi(i,k)        &  
                               + Am1(i,k)*phi(jc1,k)      &
                               + Am2(i,k)*phi(jc2,k)      &
                               + Am3(i,k)*phi(jc3,k)      &
                               + AmT(i,k)*phi(i,k+1)      &       
                               + AmB(i,k)*phi(i,k-1))     &
!                              ----------------------------
                               + Bmv1T(i,k)*phiv(jv1,k)   &
                               + Bmv2T(i,k)*phiv(jv2,k)   &
                               + Bmv3T(i,k)*phiv(jv3,k)   &
                               + Bmv1B(i,k)*phiv(jv1,k-1) &
                               + Bmv2B(i,k)*phiv(jv2,k-1) &
                               + Bmv3B(i,k)*phiv(jv3,k-1) &
!                              ----------------------------
                               + Nm(i,k) 
               enddo
            enddo

         endif !(UseDiffMethod)

      ENDIF !(TestDimension)

!*********************************************************************!
!                                                                     !
!                                 Error                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Exact solution: diffusion term             |
!     |________________________________________________________|

!     _________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do i=1,N_CELL
           do k=1,NZ 
              x = xc(i)
              y = yc(i)
              funExact(i,k) = funDiffExam2D(x,y)
              funExact(i,k) = funExact(i,k)*AreaCell(i)
           enddo
         enddo
!     _________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 
         do i=1,N_CELL
           do k=1,NZ 
              x = xc(i)
              y = yc(i)
              z = sig(k)
              if (k.eq.1) then
                 funExact(i,k) = funDiffExam3D(x,y,z)
                 funExact(i,k) = funExact(i,k)*AreaCell(i)*dsigv(1) 
              else
                 funExact(i,k) = funDiffExam3D(x,y,z)
                 funExact(i,k) = funExact(i,k)*AreaCell(i)*dsigv(k-1)
              endif
           enddo
         enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                  Approximate solution                  |
!     |________________________________________________________|

!     _________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do k=1,NZ 
            do i=1,N_CELL 
               funApprox(i,k) = diffA(i,2)
            enddo
         enddo
!     _________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 
         do k=2,NZ-1 
            do i=1,N_CELL 
               funApprox(i,k) = diffA(i,k)
            enddo
         enddo
!        -------------------------------------------
!        Solution (top & bottom)
         do i=1,N_CELL
            funApprox(i,1)  = funExact(i,1)
            funApprox(i,NZ) = funExact(i,NZ)
         enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                Calculation of the errors               |
!     |________________________________________________________|

      maxErrorA = 0.0d0
      maxErrorR = 0.0d0
      sumErrorA = 0.0d0
      sumErrorR = 0.0d0
      do k=1,NZ 
         do i=1,N_CELL0 
            ErrorA(i,k) = abs(funApprox(i,k)-funExact(i,k))
            ErrorR(i,k) = abs(funExact(i,k))
            maxErrorA = max(maxErrorA,ErrorA(i,k))
            maxErrorR = max(maxErrorR,ErrorR(i,k))
            sumErrorA = sumErrorA+ErrorA(i,k)**2
            sumErrorR = sumErrorR+ErrorR(i,k)**2
         enddo
      enddo
      sumErrorA = dsqrt(sumErrorA/(N_CELL0*NZ))
      sumErrorR = dsqrt(sumErrorR/(N_CELL0*NZ))
      maxErrorR = maxErrorA/maxErrorR
      sumErrorR = sumErrorA/sumErrorR
!      ________________________________________________________
!     |                                                        |
!     |                      Display solution                  |
!     |________________________________________________________|

9     format(t10,a3,t14,e10.3,t26,e10.3,t38,e10.3,t50,e10.3)
8     format(t10,60a)
      write(*,8)'===================================================='
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'               TEST: Diffusion term                 '
      write(*,'(t30,a4,i3)') ' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      write(*,9) '   ',maxErrorA,sumErrorA,maxErrorR,sumErrorR
      write(*,8)'----------------------------------------------------'
      write(*,8)'===================================================='
      write(*,*) '                                                   '

!*********************************************************************!
!                                                                     !
!                        Save tecplot solutions                       !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                      Cell-centers                      |
!     |________________________________________________________|

      IF (TestDimension.eq.2) THEN
         do k=1,NZ 
            do i=1,N_CELL
               SolAppro(i,k) = funApprox(i,k)/AreaCell(i)
               SolExact(i,k) = funExact(i,k)/AreaCell(i)
               SolError(i,k) = ErrorA(i,k)
           enddo
         enddo
      ELSEIF (TestDimension.eq.3) THEN           
         do i=1,N_CELL
            do k=1,NZ
               if (k.eq.1) then
                  Volume = AreaCell(i)*dsigv(1)
               else
                  Volume = AreaCell(i)*dsigv(k-1)
               endif
               SolAppro(i,k) = funApprox(i,k)/Volume
               SolExact(i,k) = funExact(i,k)/Volume
               SolError(i,k) = ErrorA(i,k)
            enddo
         enddo
       ENDIF

!      ________________________________________________________
!     |                                                        |
!     |                    Vertex Approximation                |
!     |________________________________________________________|

!     __________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do i=1,N_CELL
            phi2D(i) = SolAppro(i,3)
         enddo
!        ------------------------------------------------
!        Interpolation of the inside vertex points 
         call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                              phi2D,xc,yc,No_cp,nbe)
!        ------------------------------------------------
!        Boundary Conditions of the vertex points
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
               x = xv(nv)
               y = yv(nv)
    	       phiv2D(nv) = funDiffExam2D(x,y)
	    endif 
         enddo
!        ------------------------------------------------
!        3D format     
         do k=1,NZ-1  
            do nv=1,N_VERT
               SolApprov(nv,k) = phiv2D(nv)
            enddo
         enddo
!     __________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 

         do k=1,NZ 
            do i=1,N_CELL
               SolAppro(i,k) = SolAppro(i,k)
           enddo
         enddo
!        -------------------------------------------
!        Interpolation
         call interpolation3D(SolApprov,xv,yv,sigv,dsigv,No_vp,nbev,&
                              SolAppro,xc,yc,sig,dsig,No_cp)
!        -------------------------------------------
!        Boundary conditions
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
!           __________
!           Horizontal     
            if (nbev(nv).ne.0) then
               do k=1,NZ-1
                  z = sigv(k)                                        
    	          SolApprov(nv,k) = funDiffExam3D(x,y,z)  
               enddo
	    endif
!           __________
!           Vertical 
            z = sigv(1)
            SolApprov(nv,1)    = funDiffExam3D(x,y,z) 
            z = sigv(NZ-1)  
            SolApprov(nv,NZ-1) = funDiffExam3D(x,y,z) 
         enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                  Vertex Exact solution                 |
!     |________________________________________________________|

!     _________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do nv=1,N_VERT
           do k=1,NZ-1 
              x = xv(nv)
              y = yv(nv)
              SolExactv(nv,k) = funDiffExam2D(x,y)
           enddo
         enddo
!     _________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 
         do nv=1,N_VERT
           do k=1,NZ-1 
              x = xv(nv)
              y = yv(nv)
              z = sigv(k)
              SolExactv(nv,k) = funDiffExam3D(x,y,z)
           enddo
         enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                        Vertex Error                    |
!     |________________________________________________________|

!     __________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do i=1,N_CELL
            phi2D(i) = SolError(i,3)
         enddo
!        ------------------------------------------------
!        Interpolation of the inside vertex points 
         call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                              phi2D,xc,yc,No_cp,nbe)
!        ------------------------------------------------
!        Boundary Conditions of the vertex points
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
    	       !phiv2D(nv) = 0.0d0
	    endif 
         enddo
!        ------------------------------------------------
!        3D format     
         do k=1,NZ-1  
            do nv=1,N_VERT
               SolErrorv(nv,k) = phiv2D(nv)
            enddo
         enddo
!     __________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 
!        -------------------------------------------
!        Interpolation
         call interpolation3D(SolErrorv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              SolError,xc,yc,sig,dsig,No_cp)
!        -------------------------------------------
!        Boundary conditions
         do nv=1,N_VERT
!           __________
!           Horizontal  
            if (nbev(nv).ne.0) then
               do k=1,NZ-1                                                   
    	          SolErrorv(nv,k) = 0.0d0
               enddo
	    endif
!           __________
!           Vertical 
            SolErrorv(nv,1)    = 0.0d0
            SolErrorv(nv,NZ-1) = 0.0d0
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(ax,ay,az,bx,by,bz,fx,fy,fz)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '<<<<< End   subroutine: TEST DIFFUSION'
      write(*,*) ' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF testDiffusion                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
