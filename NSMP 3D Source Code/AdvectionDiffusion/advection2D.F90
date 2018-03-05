!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE advection2D(Am0,Am1,Am2,Am3,AmG,&
                             uu,vv,              &
                             phi,xc,yc,No_cp,xv,yv,No_vp,nbe )

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the advection contribution to the   !
!    general linear system. We have two options to the final values   !
!    of Am and AmG, depending if we choose between an implicit or     !
!    explicit scheme (cppdefs.h).                                     !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  | <--> Am1   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  !
!  | <--> Am2   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  !
!  | <--> Am3   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  !
!  | <--> AmG   |(N_CELL0,NZ)| Vector with Gradient terms          |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    | Description                         |  !
!  |____________|____________|_____________________________________|  !
!  | --> uu     |(N_CELL,NZ) | Velocity component u                |  !
!  | --> vv     |(N_CELL,NZ) | Velocity component v                |  !
!  |____________|____________|_____________________________________|  !
!  | --> phi    |(N_CELL,NZ) | Function phi at the center element  |  !
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Type of boundary cell (inside or bc)|  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !
!  |_______________|_______________________________________________|  !
!  | FluxLimiter   | (N_CELL0,NZ) Flux limiter (Vidovic 2006)      |  !
!  | FluxLimiterMin| Minimum value for the flux limiter            |  !
!  | FluxLimiterpos| Positive side of the flux limiter             |  !
!  | FluxLimiterneg| Negative side of the flux limiter             |  !
!  | FluxLimiterij | Pre-final value of the flux limiter           |  !
!  | phimin        | Maximum phi of all cell-centers related       |  !
!  | phimax        | Maximum phi of all cell-centers related       |  !
!  | d1,d2         | auxiliar diference values of phi              |  !
!  |_______________|_______________________________________________|  !
!  | dphidx        |(N_CELL,NZ) d(phi)/dx   = gradient component x |  !
!  | dphidy        |(N_CELL,NZ) d(phi)/dy   = gradient component y |  !
!  | dphidsig      |(N_CELL,NZ) d(phi)/dsig = gradient component z |  !
!  |_______________|_______________________________________________|  !
!  | C01           |(N_CELL0,NZ) Mass flux of horizontal neigh. 1  |  !
!  | C02           |(N_CELL0,NZ) Mass flux of horizontal neigh. 2  |  !
!  | C03           |(N_CELL0,NZ) Mass flux of horizontal neigh. 3  |  !
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | C0j(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | C0pos,C0neg   | Positive and negative values of the mass flux |  !
!  | GF            | Gradient (rhs contribution)                   |  !
!  | GFi,GFj       | Gradient at the cell center & neigborn        |  !
!  | GFpos,GFneg   | Positive and negative values of the gradient  |  !
!  |_______________|_______________________________________________|  !
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
!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!     ____________________________________
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|

      real*8, dimension(:) :: Am0(N_CELL0)
      real*8, dimension(:) :: Am1(N_CELL0)
      real*8, dimension(:) :: Am2(N_CELL0)
      real*8, dimension(:) :: Am3(N_CELL0)
      real*8, dimension(:) :: AmG(N_CELL0)
!     --------------------------------------
      real*8, dimension(:) :: uu(N_CELL)
      real*8, dimension(:) :: vv(N_CELL)
      real*8, dimension(:) :: phi(N_CELL)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0)
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8 ,dimension(:),allocatable :: dphidx
      real*8 ,dimension(:),allocatable :: dphidy
      real*8 :: sumfx,sumfy,deter
!     --------------------------------------
      real*8 ,dimension(:,:) :: Uij(N_CELL0,3)
      real*8 :: nxLength,nyLength
      real*8 :: u0j,v0j
!     --------------------------------------
      real*8 ,dimension(:),allocatable :: FluxLimiter
      real*8 :: FluxLimiterMin
      real*8 :: FluxLimiterpos
      real*8 :: FluxLimiterneg
      real*8 :: FluxLimiterij
      real*8 :: phimin
      real*8 :: phimax
      real*8 :: d1,d2
!     --------------------------------------
      real*8 ,dimension(:) :: aa(1:3)
      real*8 :: aa0,aaB,aaT
      real*8 :: Ue,Uijneg,Uijpos
      real*8 :: GF,GFi,GFj
      real*8 :: GFpos,GFneg
      integer:: jc,elem
!     --------------------------------------
      real*8 :: x,y,dfdxExam2D,dfdyExam2D,funSolExam2D
!     ---------------------------------------
#     ifdef KeyAdvCenter
      integer:: jj,jv1,jv2
      real*8 :: xm1,xm2,xm3,xm4
      real*8 :: ym1,ym2,ym3,ym4
      real*8 :: xmo,ymo
      real*8 :: xoe,yoe,doe
      real*8 :: xie,yie,die
      real*8 :: xje,yje,dje
      real*8 :: ddx, ddy
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Advection2D'
#     ifdef KeyAdvCenter
         write(*,'(t15,60a)'), '<---- Use   subroutine: KeyAdvCenter'
#     else
         write(*,'(t15,60a)'), '<---- Use   subroutine: KeyAdvUpwind'
#     endif       
#     endif
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

      allocate(dphidx(N_CELL),       &
               dphidy(N_CELL),       &
               FluxLimiter(N_CELL))


!*********************************************************************!
!                                                                     !
!                                Gradient                             !
!                                                                     !
!*********************************************************************!
      do i=1,N_CELL
         dphidx(i)=0.0d0
         dphidx(i)=0.0d0
      enddo

      do i=1,N_CELL0
!        __________________________________________________
!        Inside & boundary cell-centers
	 sumfx = 0.0d0
	 sumfy = 0.0d0
	 do j=1,3
	    jc = No_cp(i,j)
	    sumfx = sumfx + dxCC(i,j)*(phi(jc)-phi(i))
	    sumfy = sumfy + dyCC(i,j)*(phi(jc)-phi(i))
	 enddo
	 deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
         dphidx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	 dphidy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
!        __________________________________________________
!        Outside cell-centers  (ficticious)
!        ----------------------------------
!        Wall
	 if (nbe(i).eq.1) then
	    do j=1,3
	       nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
	       if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
	       if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
	          dphidx(nc) = dphidx(i)
	          dphidy(nc) = dphidy(i)
	       endif
	    enddo
         endif
!        ----------------------------------
!        Discharge normal
	 if (nbe(i).eq.2) then
	    do j=1,3
	       nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
	       if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
	       if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
	          dphidx(nc) = dphidx(i)
	          dphidy(nc) = dphidy(i)
	       endif
	    enddo
         endif
!        ----------------------------------
!        Water level
	 if (nbe(i).eq.3) then
	    do j=1,3
	       nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
	       if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
	       if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
	          dphidx(nc) = dphidx(i)
	          dphidy(nc) = dphidy(i)
	       endif
	    enddo
	 endif
      enddo
!      ________________________________________________________
!     |                                                        |
!     |               Use exact gradient (test case)           |
!     |________________________________________________________|

!      do i=1,N_CELL
!         x = xc(i)
!         y = yc(i)
!         dphidx(i)=dfdxExam2D(x,y)
!         dphidy(i)=dfdyExam2D(x,y)
!      enddo

!*********************************************************************!
!                                                                     !
!                        Mass flux calculation                        !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL0
         do j=1,3
!           ---------------------------------
!           Index of the neighbor cell-center
	    jc = No_cp(i,j)
!           ---------------------------------
!           Velocity at the faces
            u0j = 0.5d0*(uu(i)+uu(jc))
            v0j = 0.5d0*(vv(i)+vv(jc))
!           ---------------------------------
!           Mass fluxes of the cell
            nxLength =  dyVV(i,j)
            nyLength = -dxVV(i,j)
            Uij(i,j) = u0j*nxLength+v0j*nyLength
	 enddo
      enddo

!*********************************************************************!
!                                                                     !
!                       Flux limiter calculation                      !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Using the LED technique                |
!     |________________________________________________________|

#     ifdef KeyFluxLimiter
         do i=1,N_CELL0
!           -------------------------------
!           Maximum & minimum value of phi
            phimax = phi(i)
	    phimin = phi(i)
            do j=1,3
               jc = No_cp(i,j)
	       phimax = max(phimax,phi(jc))
	       phimin = min(phimin,phi(jc))
	    enddo
!           -------------------------------
!           Flux limiter calculation
            FluxLimiterMin = 1.0d5
	    do j=1,3
               GFi   = dphidx(i)*dxCE(i,j)+dphidy(i)*dyCE(i,j)
	       GFpos = 0.5d0*(GFi+abs(GFi))
  	       GFneg = 0.5d0*(GFi-abs(GFi))
	       if(dabs(GFi).lt.1.0d-10) then
	          FluxLimiterij = 1.0d0
	       else
	          d1 = (phimax-phi(i))/GFi
		  d2 = (phimin-phi(i))/GFi
		  FluxLimiterpos = min(1.0d0,d1)
		  FluxLimiterneg = min(1.0d0,d2)
		  FluxLimiterij  =  FluxLimiterpos*GFpos &
                                  + FluxLimiterneg*GFneg
	       endif
	       FluxLimiterMin = min(FluxLimiterMin,FluxLimiterij)
	    enddo
	    FluxLimiter(i) = FluxLimiterMin
         enddo
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                    No Flux limiter                     |
!     |________________________________________________________|

#     ifndef KeyFluxLimiter
         do i=1,N_CELL
            FluxLimiter(i) = 1.0d0
         enddo
#     endif

     do i=1, N_CELL0
!*********************************************************************!
!                                                                     !
!                       Advection contributions                       !
!                                                                     !
!*********************************************************************!
#   ifndef KeyAdvCenter  
!         __________________________________________________
!        |                                                  |
!        |                  Matrix & rhs                    |
!        |__________________________________________________|

         aa0 = 0.0d0
         GF  = 0.0d0
         do j=1,3
!           ---------------------------------
!           Number of the neighbor cell-center
	    jc = No_cp(i,j)
!           ---------------------------------
!           Mass fluxes of the cell
            Ue = Uij(i,j)
	    Uijpos =  0.5d0*(Ue+abs(Ue))
	    Uijneg = -0.5d0*(Ue-abs(Ue))
!           ---------------------------------
!           Matrix coeffients
	    aa0   =  Uijpos + aa0
	    aa(j) = -Uijneg
!           ---------------------------------
!           Gradient at of the cell & the neighbor
	    GFi =  dphidx(i)*(xme(i,j)-xc(i)) &
                 + dphidy(i)*(yme(i,j)-yc(i))

	    GFj =  dphidx(jc)*(xme(i,j)-xc(jc)) &
                 + dphidy(jc)*(yme(i,j)-yc(jc))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uijpos*GFi*FluxLimiter(i)  &
                    - Uijneg*GFj*FluxLimiter(jc)
	 enddo
#    else

!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!           Second Order Central Convection Scheme
!                          Xin BAI
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            __________________________________________________
!           |                                                  |
!           |                  Matrix & rhs                    |
!           |__________________________________________________|

!           ___________________________________________________
!           Horizontal neighbors
!           ------------------------------------
!           Assignation of mass fluxes C
            aa0 = 0.0d0
            GF  = 0.0d0
            do j=1,3
               jj=j+1
	       if (jj.gt.3) jj=jj-3
!              ---------------------------------
!              Number of the neighbor cell-center
	       jc = No_cp(i,j)
	       Ue = Uij(i,j)
!              ---------------------------------
!              Number of the neighbor vertex-point
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              ---------------------------------
!              calculate the intersection point
            xm1 = xc(i)
            xm2 = xc(jc)
            xm3 = xv(jv1)
            xm4 = xv(jv2)
            ym1 = yc(i)
            ym2 = yc(jc)
            ym3 = yv(jv1)
            ym4 = yv(jv2)
            xmo = (xm1*ym2-ym1*xm2)*(xm3-xm4) &
                 -(xm1-xm2)*(xm3*ym4-ym3*xm4)
            ymo = (xm1*ym2-ym1*xm2)*(ym3-ym4) &
                 -(ym1-ym2)*(xm3*ym4-ym3*xm4)
            ddx = (xm1-xm2)*(ym3-ym4) &
                 -(ym1-ym2)* (xm3-xm4)
            xmo = xmo/ddx 
            ymo = ymo/ddx
            xoe = xme(i,j)-xmo
            yoe = yme(i,j)-ymo
            doe = dsqrt(xoe**2 + yoe**2)
            xie = xmo-xm1
            yie = ymo-ym1
            die = dsqrt(xie**2 + yie**2)
            xje = xmo-xm2
            yje = ymo-ym2
            dje = dsqrt(xje**2 + yje**2)
!              ---------------------------------
!              Mass fluxes of the cell
	       Uijpos =  Ue*dje/(die+dje)
	       Uijneg =  Ue*die/(die+dje)
!              ---------------------------------
!              Matrix coeffients
	       aa0   =  Uijpos + aa0
	       aa(j) =  Uijneg
!              ---------------------------------
!              Gradient at of the cell & the neighbor
	       GFi =  0.5d0*(dphidx(i)*xoe &
                    + dphidy(i)*yoe)

	       GFj =  0.5d0*(dphidx(jc)*xoe &
                    + dphidy(jc)*yoe)                    
!              ---------------------------------
!              Gradient coefficient (rhs)
!               GF = GF + Ue*(GFi +GFj)
                GF = GF + 0.0
	    enddo
#     endif

!         __________________________________________________
!        |                                                  |
!        |                  Contributions                   |
!        |__________________________________________________|

!        --------------------------------------
!        Matrix
         Am0(i) = Am0(i) + aa0
         Am1(i) = Am1(i) + aa(1)
         Am2(i) = Am2(i) + aa(2)
         Am3(i) = Am3(i) + aa(3)
!        --------------------------------------
!        Known part to rhs
         AmG(i) = AmG(i) + GF
      enddo

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(dphidx,dphidy, &
                 FluxLimiter)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Advection 2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF ADVECTION  2D                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
