!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 1 !        AUXILIAR: Function eta (free surface) & h (depth)        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                     Free surface: eta                  |
!     |________________________________________________________|

      function funEta(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,funEta
      real*8  :: pi,delta 
      integer :: i
!     ----------------------------------------
!     Function
      pi = 4.0d0*datan(1.d0)
      delta = 1.0d0
!     ----------------------------------------
      !funEta  = delta*dsin(pi*(x+1)/2)*dsin(pi*(y+1)/2) + 1.0d0
      funEta  = 0.1d0*dsin(2*pi*x)*dsin(4*pi*y)+0.25d0
      !funEta  = 0.0d0

      return
      end function funEta
!     __________________________________________________________

      function funFreeSurface(x,y,t)    

      implicit none
      real*8 :: x,y,t,funFreeSurface
      real*8 :: pi,L,D,A,k,TT,s
!     ----------------------------------------

      pi = 4.0d0*datan(1.d0)
      L  = 20.0d0  ! Basin length
      D  = 10.0d0  ! Basin depth
      A  = 0.1d0   ! Amplitude of the standing wave
      k  = pi/D
      TT = 3.59d0  ! Wave period approximation
      s  = 2*pi/TT ! s^2 = g*k*tanh(k*D)

      funFreeSurface  = A*dcos(k*pi)*dcos(s*t)

      return
      end function funFreeSurface

!      ________________________________________________________
!     |                                                        |
!     |                        Depth: h                        |
!     |________________________________________________________|

      function funh(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,funh
      integer :: i
!     ----------------------------------------
!     Function   
      funh = -0.2d0*(sin(3.14*(x+y))-4.0d0)
      !funh =  1.0d0

      return
      end function funh

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 2 !           AUXILIAR: Functions of 2D examples                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                         Exam2D                         |
!     |________________________________________________________|

      function funSolExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funSolExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funSolExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funSolExam2D = (x**2-1)*(y**2-1)
      elseif (FunctionExample.eq.2) then
          funSolExam2D = dsin(pi*x)*dsin(pi*y)
      endif 

      return
      end function funSolExam2D
!      ________________________________________________________
!     |                                                        |
!     |                     Gradient Exam2D                    |
!     |________________________________________________________|
!     ________________________________________________________
!     dfdx 

      function dfdxExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,dfdxExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdxExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdxExam2D = (2.0d0*x)*(y**2-1)
      elseif (FunctionExample.eq.2) then
          dfdxExam2D = pi*dcos(pi*x)*dsin(pi*y) 
      endif 

      return
      end function dfdxExam2D
!     ________________________________________________________
!     dfdy

      function dfdyExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,dfdyExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.1) then
          dfdyExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdyExam2D = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdyExam2D = pi*dsin(pi*x)*dcos(pi*y)  
      endif 

      return
      end function dfdyExam2D
!      ________________________________________________________
!     |                                                        |
!     |         Neumman boundary condition of Exam2D           |
!     |________________________________________________________|

      function Neumanndfdn2D(x,y,nnx,nny)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,nnx,nny,Neumanndfdn2D
      real*8 :: dfdx,dfdy
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)
          dfdy = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y) 
          dfdy = pi*dsin(pi*x)*dcos(pi*y)  
      endif 
      Neumanndfdn2D = nnx*dfdx+nny*dfdy

      return
      end function Neumanndfdn2D
!      ________________________________________________________
!     |                                                        |
!     |               Advection term of Exam2D                 |
!     |________________________________________________________|

      function funAdvExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funAdvExam2D
      real*8 :: dfdx,dfdy
!     ----------------------------------------
!     Function
      if  (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)
          dfdy = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y) 
          dfdy = pi*dsin(pi*x)*dcos(pi*y)  
      endif 
      funAdvExam2D = dfdx+dfdy

      return
      end function funAdvExam2D
!      ________________________________________________________
!     |                                                        |
!     |               Diffusion term of Exam2D                 |
!     |________________________________________________________|

      function funDiffExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funDiffExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funDiffExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funDiffExam2D = 2.0d0*((y*y-1)+(x*x-1))
      elseif (FunctionExample.eq.2) then
          funDiffExam2D = -2.0d0*pi*pi*dsin(pi*x)*dsin(pi*y) 
      endif 

      return
      end function funDiffExam2D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 4 !              AUXILIAR: Functions of 3D examples                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                         Exam3D                         |
!     |________________________________________________________|

      function funSolExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funSolExam3D,c
      integer, parameter :: n = 1
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funSolExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funSolExam3D = -(x**2-1)*(y**2-1)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          !funSolExam3D = dsin(pi*x)*dsin(pi*y)*dsin(pi*z)
          !funSolExam3D = dsin(pi*(x+1)/2)*dsin(pi*(y+1)/2)*dsin(pi*(z+1)/2)
          funSolExam3D = dsin(n*pi*x)*dsin(n*pi*y)*dsin(n*pi*z)/(n*n)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          funSolExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function funSolExam3D

!      ________________________________________________________
!     |                                                        |
!     |                Diffusion term of Exam3D                |
!     |           Poisson right-hand side examples             |
!     |________________________________________________________|

      function funDiffExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funDiffExam3D
      integer, parameter :: n = 1
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funDiffExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funDiffExam3D = -2.0d0*( (y**2-1)*(z**2-1) &
                                  +(x**2-1)*(z**2-1) &
                                  +(x**2-1)*(y**2-1) ) 
      elseif (FunctionExample.eq.2) then
          !funDiffExam3D = -3.0d0*pi*pi*dsin(pi*x)*dsin(pi*y)*dsin(pi*z)
          !funDiffExam3D = -3.0d0*pi*pi/4*dsin(pi*(x+1)/2)*dsin(pi*(y+1)/2)*dsin(pi*(z+1)/2)
          funDiffExam3D = -3.0d0*pi*pi*dsin(n*pi*x)*dsin(n*pi*y)*dsin(n*pi*z)
      elseif (FunctionExample.eq.3) then
          funDiffExam3D = 0.0d0
      endif 

      return
      end function funDiffExam3D

!      ________________________________________________________
!     |                                                        |
!     |                   Gradient of Exam3D                   |
!     |________________________________________________________|
!     ________________________________________________________
!     dfdx                 
      function dfdxExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdxExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdxExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdxExam3D = (2.0d0*x)*(y**2-1)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          dfdxExam3D = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z) 
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdxExam3D = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdxExam3D
!     ________________________________________________________
!     dfdy

      function dfdyExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdyExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdyExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdyExam3D = (x**2-1)*(2.0d0*y)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          dfdyExam3D = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)  
      elseif (FunctionExample.eq.3) then
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdyExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdyExam3D
!     ________________________________________________________
!     dfdz

      function dfdzExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdzExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdzExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdzExam3D = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdzExam3D = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)  
      elseif (FunctionExample.eq.3) then
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdzExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdzExam3D
!      ________________________________________________________
!     |                                                        |
!     |          Neumman boundary condition of Exam3D          |
!     |________________________________________________________|

      function Neumanndfdn3D(x,y,z,nnx,nny,nnz)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,nnx,nny,nnz,Neumanndfdn3D,c
      real*8 :: dfdx,dfdy,dfdz
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)*(z**2-1)
          dfdy = (x**2-1)*(2.0d0*y)*(z**2-1)
          dfdz = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z)
          dfdy = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)
          dfdz = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdx = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdy = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdz = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      elseif (FunctionExample.eq.5) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      endif 
      Neumanndfdn3D = nnx*dfdx+nny*dfdy+nnz*dfdz

      return
      end function Neumanndfdn3D
!      ________________________________________________________
!     |                                                        |
!     |                Advection term of Exam3D                |
!     |________________________________________________________|

      function funAdvExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funAdvExam3D,c
      real*8 :: dfdx,dfdy,dfdz
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)*(z**2-1)
          dfdy = (x**2-1)*(2.0d0*y)*(z**2-1)
          dfdz = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z)
          dfdy = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)
          dfdz = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdx = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdy = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdz = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 
      funAdvExam3D = dfdx+dfdy+dfdz

      return
      end function funAdvExam3D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 5 !               AUXILIAR: Time problems solutions                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                           2D                           |
!     |________________________________________________________|

      function TimeExample2D(x,y,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,t,TimeExample2D
      real*8 :: x0,y0,uu,vv,rx,ry,ra,c1,c2,L0,Gam,cx,cy
!     __________________________________________
!     Example: RK2
      if (RUNtestRK2approx.eq.1) then
          TimeExample2D = dsin(t)
      endif
!     __________________________________________
!     Example: Pure advection 2D
      if (RUNtestAdvEqn.eq.1) then
         x0 = -0.5d0*dcos(2.0d0*pi*t)
         y0 = -0.5d0*dsin(2.0d0*pi*t)
         rx = x-x0
         ry = y-y0
         ra = dsqrt(rx**2+ry**2) 
         if (ra.le.0.25d0) then 
            TimeExample2D = dcos(2.0d0*pi*ra)**2
         else
            TimeExample2D = 0.0d0
         endif
      endif
!     __________________________________________
!     Example: Pure diffusion 2D
      if (RUNtestDiffEqn.eq.1) then
         Gam =  5.0d-2 
         L0  =  7.5d-2
         c2  = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1  =  1.0d0/(4.0d0*pi*(Gam*t+L0**2))         
         rx = x-0.0d0
         ry = y-0.0d0
         TimeExample2D = c1*dexp(c2*(rx**2+ry**2)) 
      endif
!     __________________________________________
!     Example: advection-diffusion 2D
      if (RUNtestAdvDiffEqn.eq.1) then
         Gam =  3.0d-2 
         L0  =  7.5d-2
         uu  =  5.0d0
         vv  =  5.0d0
         x0  = -0.5d0
         y0  = -0.5d0
         c2  = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1  =  1.0d0/(4.0d0*pi*(Gam*t+L0**2))         
         rx = x-x0-uu*t
         ry = y-y0-vv*t
         TimeExample2D = c1*dexp(c2*(rx**2+ry**2))            
      endif

      return
      end function TimeExample2D

!      ________________________________________________________
!     |                                                        |
!     |                           3D                           |
!     |________________________________________________________|

      function TimeExample3D(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,TimeExample3D
      real*8 :: x0,y0,z0,uu,vv,ww,rx,ry,rz,ra,c1,c2,L0,Gam
!     __________________________________________
!     Example: RK2
      if (RUNtestRK2approx.eq.1) then         
          TimeExample3D = dsin(t)
      endif
!     __________________________________________
!     Example: Pure advection (need to fix) 3D
      if (RUNtestAdvEqn.eq.1) then
         c2 = 0.5d0*(1.0d0/dsqrt(2.0d0))
         x0 = -c2*dcos(2.0d0*pi*t)
         y0 =  c2*dcos(2.0d0*pi*t)
         z0 = -c2*dcos(2.0d0*pi*t)
         rx =  x-x0
         ry =  y-y0
         rz =  z-z0      
         ra = dsqrt(rx**2+ry**2+rz**2)
         if (ra.le.0.25d0) then 
            TimeExample3D = dcos(2.0d0*pi*ra)**2
         else
            TimeExample3D = 0.0d0
         endif
      endif
!     __________________________________________
!     Example: Pure diffusion 3D
      if (RUNtestDiffEqn.eq.1) then
         Gam=  5.0d-2 
         L0 =  7.5d-2
         c2 = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1 =  1.0d0/dsqrt((4.0d0*pi*(Gam*t+L0**2))**3)
         rx =  x-0.0d0
         ry =  y-0.0d0
         rz =  z-0.0d0                     
         TimeExample3D = c1*dexp(c2*(rx**2+ry**2+rz**2))  
      endif
!     __________________________________________
!     Example: advection-diffusion 3D
      if (RUNtestAdvDiffEqn.eq.4) then
         Gam=  3.0d-2 
         L0 =  7.5d-2
         uu =  5.0d0
         vv =  5.0d0
         ww =  5.0d0
         x0 = -0.5d0
         y0 = -0.5d0
         z0 = -0.5d0
         c2 = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1 =  1.0d0/dsqrt((4.0d0*pi*(Gam*t+L0**2))**3)
         rx =  x-x0-uu*t
         ry =  y-y0-vv*t
         rz =  z-z0-ww*t                     
         TimeExample3D = c1*dexp(c2*(rx**2+ry**2+rz**2))  
      endif

      return
      end function TimeExample3D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 6 !              AUXILIAR: Solution N-S equations                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                Exact velocity & pressure               |
!     |________________________________________________________|

!     ________________________________________________________
!     Velocity component u

      function funExamNSu(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSu
      real*8 :: c 
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSu =-0.5d0*c*(dsqrt(3.0d0)*dcos(x)*dsin(y)*dsin(z)&
                                       + dsin(x)*dcos(y)*dcos(z))
      return
      end function funExamNSu
!     ________________________________________________________
!     Velocity component v

      function funExamNSv(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSv
      real*8 :: c 
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSv = 0.5d0*c*(dsqrt(3.0d0)*dsin(x)*dcos(y)*dsin(z)&
                                       - dcos(x)*dsin(y)*dcos(z))
      return
      end function funExamNSv
!     ________________________________________________________
!     Velocity component w

      function funExamNSw(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSw
      real*8 :: c 
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSw = c*dcos(x)*dcos(y)*dsin(z)

      return
      end function funExamNSw
!     ________________________________________________________
!     Pressure p

      function funExamNSp(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSp
      real*8 :: uu,vv,ww,c,d 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funExamNSp = -0.5d0*(uu**2+vv**2+ww**2)

      return
      end function funExamNSp

!      ________________________________________________________
!     |                                                        |
!     |        Neumman boundary condition of NS equation "p"   |
!     |________________________________________________________|

      function NeumanndpdnNS(x,y,z,nnx,nny,nnz,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,nnx,nny,nnz,NeumanndpdnNS
      real*8 :: dfdx,dfdy,dfdz
!     ----------------------------------------
!     Function
      dfdx = 0.0d0
      dfdy = 0.0d0
      dfdz = 0.0d0
      NeumanndpdnNS = nnx*dfdx+nny*dfdy+nnz*dfdz

      return
      end function NeumanndpdnNS
!      ________________________________________________________
!     |                                                        |
!     |        Function Dp = pf^(n+1)-pf^(n) (pressure)        |
!     |________________________________________________________|

      function funExamNSDp(x,y,z,tt)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,tt,funExamNSDp
      real*8 :: uu,vv,ww,c,d,funtnp,funtn 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
!     ----------------------------------------
!     Function at time (n+1) = tt
      c = dexp(-3.0d0*tt/Re)
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtnp = -0.5d0*(uu**2+vv**2+ww**2)
!     ----------------------------------------
!     Function at time (n) = tt - dt
      c = dexp(-3.0d0*(tt-dt)/Re)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtn = -0.5d0*(uu**2+vv**2+ww**2)
!     ----------------------------------------
!     Final function
      funExamNSDp = funtnp-funtn

      return
      end function funExamNSDp
!      ________________________________________________________
!     |                                                        |
!     |     Exact right-hand side of the Poisson problem       |
!     |________________________________________________________|

      function funExamNSrhsp(x,y,z,uf,vf,wf,tt)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,tt,uf,vf,wf,funExamNSrhsp
      real*8 :: c,d 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
      real*8 :: dudx,dudy,dudz
      real*8 :: dvdx,dvdy,dvdz
      real*8 :: dwdx,dwdy,dwdz
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*tt/Re)
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)

      dudx =-0.5d0*c*(-d*sinx*siny*sinz + cosx*cosy*cosz)
      dudy =-0.5d0*c*( d*cosx*cosy*sinz - sinx*siny*cosz)
      dudz =-0.5d0*c*( d*cosx*siny*cosz - sinx*cosy*sinz)
      dvdx = 0.5d0*c*( d*cosx*cosy*sinz + sinx*siny*cosz)
      dvdy = 0.5d0*c*(-d*sinx*siny*sinz - cosx*cosy*cosz)
      dvdz = 0.5d0*c*( d*sinx*cosy*cosz + cosx*siny*sinz)
      dwdx = -c*sinx*cosy*sinz
      dwdy = -c*cosx*siny*sinz
      dwdz =  c*cosx*cosy*cosz

      funExamNSrhsp  = 3.0d0*(uf**2+vf**2+wf**2)  &
                       -( dudx**2+dudy**2+dudz**2 &
                         +dvdx**2+dvdy**2+dvdz**2 &
                         +dwdx**2+dwdy**2+dwdz**2 ) 

      return
      end function funExamNSrhsp


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        End of auxiliar functions                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
