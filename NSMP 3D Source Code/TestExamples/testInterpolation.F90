!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          TEST INTERPOLATION                         !
!                              Feb 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testInterpolation(SolApprov,SolExactv,SolErrorv, &
                                   xc,yc,sig,dsig,No_cp,nbe,      &
                                   xv,yv,sigv,dsigv,No_vp,nbev)    
   
!---------------------------------------------------------------------!
!                                                                     !
!    This program test the interpolation from the cell-centers to     !
!    the vertices of the elements in 2D and 3D                        !  
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | <--SolApprov|(N_VERT,NZ)| Approximation at the vertices       |  !
!  | <--SolExactv|(N_VERT,NZ)| Exact solution at the vertices      |  !
!  | <--SolErrorv|(N_VERT,NZ)| Errors at the vertices              |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of vertex points           |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Local varaibles:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !  
!  | phi2D       |(N_CELL0)   | Function at the cell-center 2D     |  !
!  | phiv2D      |(N_VERT)    | Approx function at the vertices 2D |  !
!  | phivExact2D |(N_VERT)    | Exact  function at the vertices 2D |  !
!  | phi         |(N_CELL0,NZ)| Function at the cell-center        |  !
!  | phiv        |(N_VERT,NZ) | Approx function at the vertices    |  !
!  | phivExact   |(N_VERT,NZ) | Exact  function at the vertices    |  !
!  | ErrorA      |(N_CELL0,NZ)| Absolute error at each point       |  !
!  | ErrorR      |(N_CELLO,NZ)| Relative error at each point       |  !
!  | MaxErrorA   | (1:5,1:3)  | Maximum absolute error             |  !
!  | MaxErrorR   | (1:5,1:3)  | Maximum realtive error             |  !
!  | DisplayThis | integer    | Tag to display the solution        |  !
!  |_____________|____________|____________________________________|  !
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

      real*8,dimension(:,:) :: SolApprov(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolExactv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolErrorv(N_VERT,NZ-1)
!     -----------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -----------------------------------
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

      real*8,dimension(:),allocatable  :: phi2D
      real*8,dimension(:),allocatable  :: phiv2D
      real*8,dimension(:),allocatable  :: phivExact2D
!     -----------------------------------
      real*8,dimension(:,:),allocatable:: phi
      real*8,dimension(:,:),allocatable:: phiv
      real*8,dimension(:,:),allocatable:: phivExact
!     -----------------------------------
      integer:: jv1,jv2,jv3,s,jc
      real*8 :: x1,x2,x3,y1,y2,y3
      real*8 :: a11,a21,a12,a22,a13,a23
      real*8 :: xmy1,xmy2,xmy3,xpy1,xpy2,xpy3
      real*8 :: cp,b1,b2
      real*8 :: z1,z2,summ
      real*8 :: xx,yy,zz 
      real*8 :: x,y,z
      real*8 :: areaT
!     -----------------------------------
      real*8,dimension(:,:),allocatable:: ErrorA
      real*8,dimension(:,:),allocatable:: ErrorR
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
      integer:: ii,jj,elem
!     -----------------------------------
      integer,parameter :: UseExactAverage = 0
      integer,parameter :: DisplayThis2D   = 1
      integer,parameter :: DisplayThis3D   = 0
!     -----------------------------------
      real*8 :: funSolExam2D,funSolExam3D
!     -----------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: fB,dfBdn,nnx,nny,nnz
      real*8 :: Neumanndfdn2D,Neumanndfdn3D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: TEST interpolation'
      print*,' '
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

      allocate(phi2D(N_CELL),&
               phiv2D(N_VERT),&
               phivExact2D(N_VERT))

      allocate(phi(N_CELL,NZ),&
               phiv(N_VERT,NZ-1),&
               phivExact(N_VERT,NZ-1))

      allocate(ErrorA(N_VERT,NZ-1),&
               ErrorR(N_VERT,NZ-1))

!*********************************************************************!
!                                                                     !
!                          Interpolation 2D                           !
!                                                                     !
!*********************************************************************!

      IF (DisplayThis2D.eq.1) THEN

!         _____________________________________________________
!        |                                                     |
!        |             Function at cell-centers                |
!        |_____________________________________________________|

!        _____________________________________________________
!        Direct
         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
            phi2D(i) = funSolExam2D(x,y)
         enddo
!        _____________________________________________________
!        Average: fun = int(phi2D)/m(A)
         IF (UseExactAverage.eq.1) THEN
             DO i=1,N_CELL0
                summ = 0.0d0
                jv1 = No_vp(i,1)
                jv2 = No_vp(i,2)
                jv3 = No_vp(i,3)
                x1 = xv(jv1)
                y1 = yv(jv1) 
                x2 = xv(jv2)
                y2 = yv(jv2) 
                x3 = xv(jv3)
                y3 = yv(jv3) 
                a11 = (x2+x3-2.0d0*x1)/2.0d0
                a21 = (y2+y3-2.0d0*y1)/2.0d0
                a12 = (x3-x2)/2.0d0
                a22 = (y3-y2)/2.0d0
                a13 = x1
                a23 = y1
                do s=1,2
                   if (s.eq.1) then        
                      xmy1 = pi*(a11-a21)        
                      xmy2 = pi*(a12-a22)        
                      xmy3 = pi*(a13-a23)
                      cp   = 1.0d0
                   elseif (s.eq.2) then
                      xmy1 = pi*(a11+a21)
                      xmy2 = pi*(a12+a22)
                      xmy3 = pi*(a13+a23)
                      cp   = -1.0d0
                   endif
                   if (abs(xmy2).le.1d-9) then
                      if (abs(xmy1).le.1d-9) then
                         summ = cp*0.5d0*dcos(xmy3)
                      else
                         b1 = 1.0d0/xmy1
                         b2 = b1*b1
                         summ = summ+cp*b1*(dsin(xmy1+xmy3)) &
                                    +cp*b2*(dcos(xmy1+xmy3)-dcos(xmy3))
                      endif
                   else
                      b1 = 1.0d0/xmy2
                      if (abs(xmy1+xmy2).le.1d-9) then
                         summ = summ + cp*0.5d0*b1*dsin(xmy3)
                      else
                         b2 = -1.0d0/(xmy1+xmy2)
                         summ = summ &
                              +cp*0.5d0*b1*b2*(dcos(xmy1+xmy2+xmy3)&
                                              -dcos(xmy3))
                      endif
                      if (abs(xmy1-xmy2).le.1d-9) then
                         summ = summ - cp*0.5d0*b1*dsin(xmy3)
                      else
                         b2 = -1.0d0/(xmy1-xmy2)        
                         summ = summ &
                              -cp*0.5d0*b1*b2*(dcos(xmy1-xmy2+xmy3)&
                                              -dcos(xmy3))
                      endif
                   endif
                enddo
                phi2D(i) = summ                
             ENDDO
         ENDIF
!         _____________________________________________________
!        |                                                     |
!        |   Boundary conditions of the cell-centers points    |
!        |_____________________________________________________|

!        _____________________________________________________
!        Dirichlet
         IF (ChooseBoundary.eq.1) THEN
            do i=1,N_CELL0
	       if (nbe(i).ne.0) then	
   	          do j=1,3
	             nc=No_cp(i,j)
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
                        fB = funSolExam2D(x,y)   
                        phi2D(nc) = 2.0d0*fB-phi2D(i)               
                     endif 
	          enddo
               endif
            enddo
!        _____________________________________________________
!        Neumann
         ELSEIF (ChooseBoundary.eq.2) THEN
            do i=1,N_CELL0
	       if (nbe(i).ne.0) then	
   	          do j=1,3
	             nc=No_cp(i,j)
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
                         nnx = normxc(i,j)
                         nny = normyc(i,j)
                         dfBdn = Neumanndfdn2D(x,y,nnx,nny)
                         phi2D(nc) = 2.0d0*dlCE(i,j)*dfBdn + phi2D(i)
                     endif 
	          enddo
               endif
            enddo
          ENDIF
!         _____________________________________________________
!        |                                                     |
!        |                    Interpolation                    |
!        |_____________________________________________________|
        
         call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                              phi2D,xc,yc,No_cp,nbe)

!         _____________________________________________________
!        |                                                     |
!        |      Boundary condition at the vertex points        |
!        |_____________________________________________________|

         call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                         phi2D,xc,yc,No_cp,nbe)

!         _____________________________________________________
!        |                                                     |
!        |             Exact solution at the vertices          |
!        |_____________________________________________________|

         do nv=1,N_VERT  
            x = xv(nv)
            y = yv(nv)
            phivExact2D(nv) = funSolExam2D(x,y) 
         enddo
!         _____________________________________________________
!        |                                                     |
!        |                        Error                        |
!        |_____________________________________________________|

!        _____________________________________________________
!        Error (Maximum norm: p=infty)         
         maxErrorA = 0.0d0
         maxErrorR = 0.0d0
         do k=1,NZ-1
            do nv=1,N_VERT  
               ErrorA(nv,k) = abs(phivExact2D(nv)-phiv2D(nv))
               ErrorR(nv,k) = abs(phivExact2D(nv)) 
               maxErrorA = max(maxErrorA,ErrorA(nv,k))
               maxErrorR = max(maxErrorR,ErrorR(nv,k))
            enddo
         enddo
         maxErrorR = maxErrorA/maxErrorR 
!        _____________________________________________________
!        Error (Frobenius norm: p=2)         
         sumErrorA = 0.0d0
         sumErrorR = 0.0d0
         k=2
         do nv=1,N_VERT  
            sumErrorA = sumErrorA +(ErrorA(nv,k))**2
            sumErrorR = sumErrorR +(ErrorR(nv,k))**2
         enddo
         sumErrorA = dsqrt(sumErrorA/N_VERT)
         sumErrorR = dsqrt(sumErrorR/N_VERT)
         sumErrorR = sumErrorA/sumErrorR 
!        _____________________________________________________
!        Save solutions
         do k=1,NZ-1 
            do nv=1,N_VERT
               SolApprov(nv,k) = phiv2D(nv) 
               SolExactv(nv,k) = phivExact2D(nv)
               SolErrorv(nv,k) = ErrorA(nv,k)
            enddo
         enddo
      ENDIF

9     format(t10,a3,t14,e10.3,t26,e10.3,t38,e10.3,t50,e10.3)
8     format(t10,60a)
      write(*,8)'===================================================='
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'             TEST: interpolation vertex             '
      write(*,'(t30,a4,i3)') ' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      if (DisplayThis2D.eq.1) then
          write(*,9) '2D:',maxErrorA,sumErrorA,maxErrorR,sumErrorR
      endif

!*********************************************************************!
!                                                                     !
!                          Interpolation 3D                           !
!                                                                     !
!*********************************************************************!


      IF (DisplayThis3D.eq.1) THEN
!         _____________________________________________________
!        |                                                     |
!        |             Function at cell-centers                |
!        |_____________________________________________________|

         do i=1,N_CELL0
            do k=2,NZ-1 
               x = xc(i)
               y = yc(i)
               z = sig(k)
               phi(i,k) = funSolExam3D(x,y,z) 
            enddo
         enddo
!         _____________________________________________________
!        |                                                     |
!        |   Boundary conditions of the cell-centers points    |
!        |_____________________________________________________|

!        _____________________________________________________
!        Dirichlet
         IF (ChooseBoundary.eq.1) THEN
!           --------------------------
!           Vertical
            do i=1,N_CELL0
               x = xc(i)
               y = yc(i)
!              ______                    
!              Bottom
               k = 1
               z = sigv(k)
               fB = funSolExam3D(x,y,z) 
               phi(i,k)  = 2.0d0*fB-phi(i,k+1)
!              ______                   
!              Top
               k = NZ
               z = sigv(k-1)
               fB = funSolExam3D(x,y,z)
               phi(i,k)= 2.0d0*fB-phi(i,k-1)
            enddo
!           --------------------------
!           Horizontal
            do i=1,N_CELL0
	       if (nbe(i).ne.0) then	
   	          do j=1,3
	            nc=No_cp(i,j)
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
                       do k=1,NZ 
                          x = xe(i,j)
                          y = ye(i,j)
                          z = sig(k)
                          fB = funSolExam3D(x,y,z)   
                          phi(nc,k) = 2.0d0*fB-phi(i,k)                    
                       enddo
                    endif 
	         enddo
               endif
            enddo
!        _____________________________________________________
!        Neumann
         ELSEIF (ChooseBoundary.eq.2) THEN
!           --------------------------
!           Vertical
            do i=1,N_CELL0
               x = xc(i)
               y = yc(i)
!              ______                    
!              Bottom
               k = 1
               z = sigv(k)
               nnx =  0.0d0
               nny =  0.0d0
               nnz = -1.0d0
               dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
               phi(i,k) = dsig(k)*dfBdn + phi(i,k+1) 
!              ______                   
!              Top
               k = NZ
               z = sigv(k-1)
               nnx = 0.0d0
               nny = 0.0d0
               nnz = 1.0d0
	       dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
               phi(i,k) = dsig(k-1)*dfBdn + phi(i,k-1)
            enddo
!           --------------------------
!           Horizontal
            do i=1,N_CELL0
	       if (nbe(i).ne.0) then	
   	          do j=1,3
	            nc=No_cp(i,j)
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
                       nnx = normxc(i,j)
                       nny = normyc(i,j)
                       nnz = 0.0d0
                       do k=1,NZ 
                          x = xe(i,j)
                          y = ye(i,j)
                          z = sig(k)
                          dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
                          phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn + phi(i,k)
                       enddo
                    endif 
	         enddo
               endif
            enddo
         ENDIF
!         _____________________________________________________
!        |                                                     |
!        |                    Interpolation                    |
!        |_____________________________________________________|

         call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                              phi,xc,yc,sig,dsig,No_cp,nbe)
!         _____________________________________________________
!        |                                                     |
!        |      Boundary condition at the vertex points        |
!        |_____________________________________________________|

         call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                         phi,xc,yc,sig,dsig,No_cp,nbe)
!         _____________________________________________________
!        |                                                     |
!        |             Exact solution at the vertices          |
!        |_____________________________________________________|

         do k=1,NZ-1 
            do nv=1,N_VERT  
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)
               phivExact(nv,k) = funSolExam3D(x,y,z)
            enddo
         enddo
!         _____________________________________________________
!        |                                                     |
!        |                        Error                        |
!        |_____________________________________________________|

!        _____________________________________________________
!        Error (Maximum norm: p=infty) 
         maxErrorA = 0.0d0
         maxErrorR = 0.0d0
         do k=1,NZ-1 
            do nv=1,N_VERT  
               ErrorA(nv,k) = abs(phivExact(nv,k)-phiv(nv,k))
               ErrorR(nv,k) = abs(phivExact(nv,k)) 
               maxErrorA = max(maxErrorA,ErrorA(nv,k))
               maxErrorR = max(maxErrorR,ErrorR(nv,k))
            enddo
         enddo
         maxErrorR = maxErrorA/maxErrorR 
!        _____________________________________________________
!        Error (Frobenius norm: p=2)         
         sumErrorA = 0.0d0
         sumErrorR = 0.0d0
         do k=1,NZ-1
            do nv=1,N_VERT  
               sumErrorA = sumErrorA +(ErrorA(nv,k))**2
               sumErrorR = sumErrorR +(ErrorR(nv,k))**2
            enddo
         enddo
         sumErrorA = dsqrt(sumErrorA/(N_VERT*(NZ-1)))
         sumErrorR = dsqrt(sumErrorR/(N_VERT*(NZ-1)))
         sumErrorR = sumErrorA/sumErrorR
!        _____________________________________________________
!        Save solutions
         do k=1,NZ-1 
           do nv=1,N_VERT
              SolApprov(nv,k) = phiv(nv,k) 
              SolExactv(nv,k) = phivExact(nv,k)
              SolErrorv(nv,k) = ErrorA(nv,k)
           enddo
         enddo
      ENDIF

      if (DisplayThis3D.eq.1) then
          write(*,9) '3D:',maxErrorA,sumErrorA,maxErrorR,sumErrorR
      endif
      write(*,8)'----------------------------------------------------'
      write(*,8)'===================================================='
      write(*,*) ''

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(phi2D,phiv2D,phivExact2D)
      deallocate(phi,phiv,phivExact)
      deallocate(ErrorA,ErrorR)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '<<<<< End   subroutine: TEST interpolation'
      write(*,*) ' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
