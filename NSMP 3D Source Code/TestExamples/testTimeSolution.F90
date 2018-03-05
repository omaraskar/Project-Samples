!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          AUXILIAR 1: Exact Solution advection-diffusion             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeSolution(ExactSol,xc,yc,sig,ttime)

!---------------------------------------------------------------------!
!                                                                     !
!     This subroutine calculates the exact solution of the problem    !
!     at any time and point.                                          !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
!     |   Declaration of variables         |
!     |____________________________________|

      real*8,dimension(:,:) :: ExactSol(N_CELL,NZ)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8 :: ttime
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: Gam,L0,uu,vv,ww,x0,y0,z0
      real*8 :: c1,c2,rx,ry,rz

!*********************************************************************!
!                                                                     !
!                            Functions                                !
!                                                                     !
!*********************************************************************!

!     ______________________________________________
!     Example: RK
      if (TimeExample.eq.1) then

          do i=1,N_CELL
             do k=1,NZ             
                ExactSol(i,k) = dsin(ttime)
             enddo
          enddo
!     ______________________________________________
!     Example: Pure advection
      elseif (TimeExample.eq.2) then

          do i=1,N_CELL
             do k=1,NZ             
                ExactSol(i,k) = 0.0d0            
             enddo
          enddo
!     ______________________________________________
!     Example: advection-diffusion
      elseif (TimeExample.eq.3) then

         Gam =  3.0d-2 
         L0  =  7.5d-2
         uu  =  5.0d0
         vv  =  5.0d0
         ww  =  5.0d0
         x0  = -0.5d0
         y0  = -0.5d0
         z0  = -0.5d0
         c2  = -1.0d0/(4.0d0*(Gam*ttime+L0**2))
!        -----------------------------
!        2D
         if (TestDimension.eq.2) then
            c1 =  1.0d0/(4.0d0*pi*(Gam*ttime+L0**2))
            do i=1,N_CELL
               do k=1,NZ             
                  rx = xc(i)-x0-uu*ttime
                  ry = yc(i)-y0-vv*ttime
                  ExactSol(i,k) = c1*dexp(c2*(rx**2+ry**2))            
               enddo
            enddo
!        -----------------------------
!         3D
         elseif (TestDimension.eq.3) then

            c1 =  1.0d0/dsqrt((4.0d0*pi*(Gam*ttime+L0**2))**3)
            do i=1,N_CELL
               do k=1,NZ               
                  rx =  xc(i)-x0-uu*ttime
                  ry =  yc(i)-y0-vv*ttime
                  rz = sig(k)-z0-ww*ttime                      
                  ExactSol(i,k) = c1*dexp(c2*(rx**2+ry**2+rz**2))  
               enddo
            enddo
         endif
      endif

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           End of auxiliar 1                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
