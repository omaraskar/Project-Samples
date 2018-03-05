!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!        AUXILIAR: Free surface Linear standing waves solution        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                     Exact free surface                 |
!     |________________________________________________________|

      function FS_funeta(x,y,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,t,FS_funeta
      real*8 :: L,W,A,kx,ky,TT,omega
      real*8, parameter :: pi = 3.14159265359d0
!     ----------------------------------------
      L  = 10.0d0 ! m
      W  = 10.0d0 ! m
      A  = 0.1d0  ! m
      kx = pi/L   
      ky = pi/W
      TT  = 3.01d0 ! s
      omega = 2*pi/TT
!     ----------------------------------------             

!     Function
      FS_funeta = A*dcos(kx*x)*dcos(ky*y)*dcos(omega*t)
      return
      end function FS_funeta

!      ________________________________________________________
!     |                                                        |
!     |                     Exact velocity                     |
!     |________________________________________________________|

!     ________________________________________________________
!     Velocity component u

      function FS_funu(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,z,t,FS_funu
      real*8 :: g,L,W,h,A,kx,ky,k,TT,omega,c
      real*8 :: pi = 3.14159265359
      g  =  9.80665d0        
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0       
      kx = pi/L   
      ky = pi/W
      k  = 0.44d0 ! = sqrt(kx**2 +ky**2)        
      TT = 3.01d0   
      omega = 2*pi/TT
!     ----------------------------------------
!     Function
      c = (A*g*kx/omega)*dcosh(k*(h+z))/dcosh(k*h)
      FS_funu = c*dsin(kx*x)*dcos(ky*y)*dsin(omega*t)
      return
      end function FS_funu
!     ________________________________________________________
!     Velocity component v

      function FS_funv(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,z,t,FS_funv
      real*8 :: g,L,W,h,A,kx,ky,k,TT,omega,c
      real*8 :: pi = 3.14159265359
      g  = 9.80665d0       
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0       
      kx = pi/L   
      ky = pi/W
      k  = 0.44d0       
      TT = 3.01d0 
      omega = 2*pi/TT
!     ----------------------------------------
!     Function
      c = (A*g*ky/omega)*dcosh(k*(h+z))/dcosh(k*h)
      FS_funv = c*dcos(kx*x)*dsin(ky*y)*dsin(omega*t)
      return
      end function FS_funv
!     ________________________________________________________
!     Velocity component w

      function FS_funw(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,z,t,FS_funw
      real*8 :: g,L,W,h,A,kx,ky,k,TT,omega,c
      real*8 :: pi = 3.14159265359
      g  = 9.80665d0       
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0       
      kx = pi/L   
      ky = pi/W
      k  = 0.44d0      
      TT = 3.01d0 
      omega = 2*pi/TT
!     ----------------------------------------
!     Function
      c = - (A*g*k/omega)*dsinh(k*(h+z))/dcosh(k*h)
      FS_funw = c*dcos(kx*x)*dcos(ky*y)*dsin(omega*t)

      return
      end function FS_funw

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        End of auxiliar functions                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
