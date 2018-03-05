!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!         3D BOUNDARY CONDITION CELL-CENTERS WITH FREE SURFACE        !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_BCvelc(phi,xc,yc,sig,dsig,No_cp,nbe,&
                           Hpr,h,eta,tagBC)
  
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition of the velocity. The tag called "tagBC" is used to     !
!    choose between velocity components:                              !
!                          tagBC = 0 is dphi/dn=0                     !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--- phi    |(N_CELL,NZ)| Function at the cell-center         |  !   
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | ---> Hpr    |(N_CELL)   | Total depth H = eta + h             |  !
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> sig    |(NZ)       | sigma value at the cell centers     |  !
!  | ---> dsig   |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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
      real*8, dimension(:)  :: eta(N_CELL)
!     --------------------------------------
      integer :: tagBC,elem
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB
      real*8 :: z0,z1,z2,z3,f1,f2,f3,a1,a2,a3
      real*8 :: FS_funu,FS_funv,FS_funw

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: FS_BCvelc'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                             Initialization                          !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                       d(phi)/dn = 0                    |
!     |________________________________________________________|

      IF (tagBC.eq.0) THEN
         DO i=1,N_CELL0
!           __________________________________________________
!           Vertical
!           ______                    
!           Bottom
            phi(i,1) = phi(i,2)
!           ______                   
!           Top
            phi(i,NZ)= phi(i,NZ-1)
!           __________________________________________________
!           Horizontal
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
                        phi(nc,k) = phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component u             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.1) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z  = 0.5d0*(sig(1)+sig(2))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funu(x,y,z0,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funu(x,y,z0,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
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
                        z = sig(k)*Hpr(i)-h(i)
                        fB = Hpr(i)*FS_funu(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component v             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.2) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = 0.5d0*(sig(1)+sig(2))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funv(x,y,z0,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funv(x,y,z0,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
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
                        z = sig(k)*Hpr(i)-h(i)
                        fB = Hpr(i)*FS_funv(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                     enddo
                  endif 
	       enddo
            endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component w             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.3) THEN
         DO i=1,N_CELL0
!           --------------------------------------------------
!           VERTICAL
            x = xc(i)
            y = yc(i)
!           ______                    
!           Bottom
            k = 1
            z = 0.5d0*(sig(1)+sig(2))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funw(x,y,z0,time) 
            z1 = sig(1)
            z2 = sig(2)
            z3 = sig(3)
            f2 = phi(i,2)
            f3 = phi(i,3)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f1 = (fB-a2*f2-a3*f3)/a1
            !phi(i,1) = f1
            phi(i,1) = 2.0d0*fB-phi(i,2)
!           ______                   
!           Top
            k = NZ
            z = 0.5d0*(sig(NZ)+sig(NZ-1))
            z0 = z*Hpr(i)-h(i)
            fB = Hpr(i)*FS_funw(x,y,z0,time)
            z1 = sig(NZ-2)
            z2 = sig(NZ-1)
            z3 = sig(NZ)
            f1 = phi(i,NZ-2)
            f2 = phi(i,NZ-1)
            a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
            a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
            a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
            f3 = (fB-a1*f1-a2*f2)/a3
            !phi(i,NZ) = f3
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
                        z = sig(k)*Hpr(i)-h(i)
                        fB = Hpr(i)*FS_funw(x,y,z,time)   
                        phi(nc,k) = 2.0d0*fB-phi(i,k)                    
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
           print*, '      <----   End subroutine: FS_BCvelc'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          VELOCITY VERTEX BOUNDARY CONDITION WITH FREE SURFACE       !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_BCvelv(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                           phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition    !
!    of the velocity components. The tag called "tagBC" is used to    !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- phiv   |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | ---> Hprv   |(N_VERT)   | Total depth Hv = etav + hv          |  !
!  | ---> xv,yv  |(N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv   |(NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv  |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp  |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | ---> nbev   |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !  
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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
      real*8, dimension(:)  :: eta(N_CELL)
!     --------------------------------------
      integer :: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: ic1,ic2,jc
      real*8 :: x,y,z,fB,dfBdn
      real*8 :: FS_funu,FS_funv,FS_funw
!     --------------------------------------
      real*8,dimension(:,:):: dphidx(N_CELL,NZ)
      real*8,dimension(:,:):: dphidy(N_CELL,NZ)
      real*8 :: sumfx,sumfy,deter
!     --------------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
!     --------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: FS_BCvelv'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!              Dirichlet Boundary Condition for the velocity          !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                       d(phiv)/dn = 0                   |
!     |________________________________________________________|

      IF (tagBC.eq.0) THEN
!        _____________________________________________________
!        Vertical
         do nv=1,N_VERT
!           _________                 
!           Bottom
            k  = 1
            dfBdn = 0.0d0
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
            dfBdn = 0.0d0
            f1 = phiv(nv,k-1) 
            f2 = phiv(nv,k-2) 
            h1 = sigv(k-1)-sigv(k)
            h2 = sigv(k-2)-sigv(k)
            deno = (h1*h2*h2-h2*h1*h1)
            f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
            phiv(nv,k)= f0
         enddo
!        _____________________________________________________
!        Horizontal
!        -----------------------------------
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
!        ----------------------------------- 
!        Neumann BC approximation  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
               do k=1,NZ-1
                  dfBdn = 0.0d0
!                 ----------------
!                 Function approximations
                  ic1 = icxn1(nv)
                  ic2 = icxn2(nv)
                  f1 = phi(ic1,k) + dphidx(ic1,k)*(xn1(nv)-xc(ic1)) &
                                  + dphidy(ic1,k)*(yn1(nv)-yc(ic1))   
                  f2 = phi(ic2,k) + dphidx(ic2,k)*(xn2(nv)-xc(ic2)) &
                                  + dphidy(ic2,k)*(yn2(nv)-yc(ic2))
!                 ----------------
!                 Approx at the vertex point
                  h1 = 1.0d0*dn(nv)
                  h2 = 2.0d0*dn(nv)
                  deno = (h1*h2*h2-h2*h1*h1)
                  f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2*h2-h1*h1)
                  phiv(nv,k) = f0 
               enddo
            endif           
         enddo
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component u             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.1) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)*Hprv(nv)-hv(nv)
            phiv(nv,k) = Hprv(nv)*FS_funu(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)*Hprv(nv)-hv(nv)  
            phiv(nv,k) = Hprv(nv)*FS_funu(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  phiv(nv,k) = Hprv(nv)*FS_funu(x,y,z,time) 
               enddo
           endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component v             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.2) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)*Hprv(nv)-hv(nv)
            phiv(nv,k) = Hprv(nv)*FS_funv(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)*Hprv(nv)-hv(nv)  
            phiv(nv,k) = Hprv(nv)*FS_funv(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  phiv(nv,k) = Hprv(nv)*FS_funv(x,y,z,time) 
               enddo
           endif
         ENDDO
!      ________________________________________________________
!     |                                                        |
!     |             Dirichlet velocity component w             |
!     |________________________________________________________|

      ELSEIF (tagBC.eq.3) THEN
         DO nv=1,N_VERT
!           --------------------------------------------------
!           VERTICAL
            x = xv(nv)
            y = yv(nv)
!           _________                 
!           Bottom
            k = 1
            z  = sigv(k)*Hprv(nv)-hv(nv)
            phiv(nv,k) = Hprv(nv)*FS_funw(x,y,z,time) 
!           _________                
!           Top
            k = NZ-1
            z = sigv(k)*Hprv(nv)-hv(nv)  
            phiv(nv,k) = Hprv(nv)*FS_funw(x,y,z,time)
!           --------------------------------------------------
!           HORIZONTAL
            if (nbev(nv).ne.0) then
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  phiv(nv,k) = Hprv(nv)*FS_funw(x,y,z,time) 
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
           print*, '      <----   End subroutine: FS_BCvelv'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                                  END                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
