!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          TESTS: EDGE VALUES                         !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testEdgeValues(SolAppro,SolExact,SolError, &
                                xc,yc,sig,dsig,No_cp,nbe,   &
                                xv,yv,sigv,dsigv,No_vp,nbev)    
   
!---------------------------------------------------------------------!
!                                                                     !
!    This program test the approximation of the edges values using    !
!    the gradient approximation.                                      !    
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name      |    Size   | Description                       |  !  
!  |_______________|___________|___________________________________|  !  
!  | <-- SolApprov |(N_CELL,NZ)| Approximation at the cell-centers |  !
!  | <-- SolExactv |(N_CELL,NZ)| Exact solution at the cell-centers|  !
!  | <-- SolErrorv |(N_CELL,NZ)| Errors at the vertices            |  !
!  |_______________|___________|___________________________________|  !
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
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
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

      real*8,dimension(:,:),allocatable :: phi
      real*8,dimension(:,:),allocatable :: Averagephi
!     -------------------------------------
      real*8,dimension(:,:),allocatable  :: dfdx
      real*8,dimension(:,:),allocatable  :: dfdy
      real*8,dimension(:,:),allocatable  :: dfdsig
      real*8,dimension(:,:),allocatable  :: Exactdfdx
      real*8,dimension(:,:),allocatable  :: Exactdfdy
      real*8,dimension(:,:),allocatable  :: Exactdfdsig
!     -------------------------------------
      real*8,dimension(:,:),allocatable :: funApprox
      real*8,dimension(:,:),allocatable :: funExact
      real*8,dimension(:,:),allocatable :: ErrorA1
      real*8,dimension(:,:),allocatable :: ErrorR1
      real*8,dimension(:,:) :: MaxErrorA1(1:5)
      real*8,dimension(:,:) :: MaxErrorR1(1:5)
      real*8,dimension(:,:) :: sumErrorA(1:5)
      real*8,dimension(:,:) :: sumErrorR(1:5)
!     -------------------------------------
      real*8,dimension(:,:),allocatable :: funEdge1
      real*8,dimension(:,:),allocatable :: funEdge2
      real*8,dimension(:,:),allocatable :: funEdge3
!     -------------------------------------
      real*8 :: xx,yy,zz 
      integer:: s,Nface
!     -------------------------------------
      real*8 :: x1,x2,x3,y1,y2,y3,mx1,mx2,mx3,my1,my2,my3,z1,z2,summ
      real*8 :: a11,a21,a12,a22,a13,a23,xmy1,xmy2,xmy3,xpy1,xpy2,xpy3
      real*8 :: cp,b1,b2
      integer:: jv1,jv2,jv3,jc,elem
!     -------------------------------------
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: Neumanndfdn2D,Neumanndfdn3D
      real*8 :: x,y,z,fB,dfBdn,nnx,nny,nnz
!     -------------------------------------
      real*8 :: dfdxExam2D,dfdyExam2D
      real*8 :: dfdxExam3D,dfdyExam3D,dfdzExam3D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,*) ''
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: TEST Edge Values'
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

      allocate(phi(N_CELL,NZ),&
               Averagephi(N_CELL,NZ))

      allocate(dfdx(N_CELL,NZ),&
               dfdy(N_CELL,NZ),&
               dfdsig(N_CELL,NZ))

      allocate(Exactdfdx(N_CELL,NZ),&
               Exactdfdy(N_CELL,NZ),&
               Exactdfdsig(N_CELL,NZ))

      allocate(funApprox(N_CELL0,NZ), &
               funExact(N_CELL0,NZ),  &
               ErrorA1(N_CELL0,NZ),   &
               ErrorR1(N_CELL0,NZ))

      allocate(funEdge1(N_CELL0,NZ), &
               funEdge2(N_CELL0,NZ), &
               funEdge3(N_CELL0,NZ))


!*********************************************************************!
!                                                                     !
!                              Funtion 2D                             !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.2) THEN
!         _____________________________________________________
!        |                                                     |
!        |            2D Function at cell-centers              |
!        |_____________________________________________________|

         do i=1,N_CELL0
            do k=1,NZ 
               x = xc(i)
               y = yc(i)
               phi(i,k) = funSolExam2D(x,y)
            enddo
         enddo
!         _____________________________________________________
!        |                                                     |
!        |   Boundary conditions of the cell-centers points    |
!        |_____________________________________________________|

!        _____________________________________________________
!        Dirichlet
         if (ChooseBoundary.eq.1) then
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
                        do k=1,NZ 
                           phi(nc,k) = 2.0d0*fB-phi(i,k)
                        enddo
                     endif 
	          enddo
               endif
            enddo
!        _____________________________________________________
!        Neumann
         elseif (ChooseBoundary.eq.2) then
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
                         x = xme(i,j)
                         y = yme(i,j)
                         nnx = normxc(i,j)
                         nny = normyc(i,j)
                         dfBdn = Neumanndfdn2D(x,y,nnx,nny)
                         do k=1,NZ 
                            phi(nc,k) = 2.0d0*dlCE(i,j)*dfBdn+phi(i,k)
                         enddo
                     endif 
	          enddo
               endif
            enddo
         endif
      ENDIF

!*********************************************************************!
!                                                                     !
!                              Funtion 3D                             !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.3) THEN
!         _____________________________________________________
!        |                                                     |
!        |             Function at cell-centers                |
!        |_____________________________________________________|

         do i=1,N_CELL0
            do k=1,NZ 
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
         if (ChooseBoundary.eq.1) then
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
         elseif (ChooseBoundary.eq.2) then
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
         endif
      ENDIF

!*********************************************************************!
!                                                                     !
!                        Exact Gradient function 3D                   !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.2) THEN
         do k=1,NZ 
            do i=1,N_CELL0
               x = xc(i)
               y = yc(i)
               Exactdfdx(i,k)   = dfdxExam2D(x,y)
               Exactdfdy(i,k)   = dfdyExam2D(x,y)
               Exactdfdsig(i,k) = 0.0d0
            enddo
         enddo
      ELSEIF (TestDimension.eq.3) THEN    
         do k=1,NZ 
            do i=1,N_CELL0
               x = xc(i)
               y = yc(i)
               z = sig(k)
               Exactdfdx(i,k)   = dfdxExam3D(x,y,z) 
               Exactdfdy(i,k)   = dfdyExam3D(x,y,z) 
               Exactdfdsig(i,k) = dfdzExam3D(x,y,z) 
            enddo
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!               Average value of the sin function at the edge         !
!                                                                     !
!*********************************************************************!
     
!      ________________________________________________________
!     |                                                        |
!     |             Exact integral: int(fun)/m(A)              |
!     |________________________________________________________|

     do i=1,N_CELL0
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
                 summ = summ + cp*b1*(dsin(xmy1+xmy3)) &
                             + cp*b2*(dcos(xmy1+xmy3)-dcos(xmy3))
              endif
           else
              b1 = 1.0d0/xmy2
              if (abs(xmy1+xmy2).le.1d-9) then
                 summ = summ + cp*0.5d0*b1*dsin(xmy3)
              else
                 b2 = -1.0d0/(xmy1+xmy2)
                 summ = summ &
                      +cp*0.5d0*b1*b2*(dcos(xmy1+xmy2+xmy3)-dcos(xmy3))
              endif

              if (abs(xmy1-xmy2).le.1d-9) then
                 summ = summ - cp*0.5d0*b1*dsin(xmy3)
              else
                 b2 = -1.0d0/(xmy1-xmy2)
                 summ = summ &
                      -cp*0.5d0*b1*b2*(dcos(xmy1-xmy2+xmy3)-dcos(xmy3))
              endif
           endif
        enddo
        
        if (TestDimension.eq.2) then
           do k=1,NZ
              Averagephi(i,k) = summ
           enddo
        elseif (TestDimension.eq.3) then
           do k=1,NZ
              z1 = -dcos(pi*sigv(k))/pi
              z2 = -dcos(pi*sigv(k+1))/pi
              Averagephi(i,k) = summ*(z2-z1)/dsigv(k)
           enddo
        endif
     enddo

      do i=1,N_CELL0
	if (nbe(i).ne.0) then	
	   do j=1,3
	      nc=No_cp(i,j)
              if (nc.gt.N_CELL0) then                 
                 do k=1,NZ                      
                    Averagephi(nc,k) = -Averagephi(i,k)                     
                 enddo
              endif 
	   enddo
         endif
      enddo

!*********************************************************************!
!                                                                     !
!               Approximate a function at each edge point             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |           Approximation Gradient function 3D           |
!     |________________________________________________________|

      call grandientLSM(dfdx,dfdy,dfdsig,&
                        phi,No_cp,nbe,sig) 

!      ________________________________________________________
!     |                                                        |
!     |     Approximation and exact function at the edge       |
!     |________________________________________________________|


      if (TestDimension.eq.2) Nface = 3    
      if (TestDimension.eq.3) Nface = 5
      DO s=1,Nface
         do k=1,NZ
            do i=1,N_CELL0
               if (s.eq.1) then
                  xx = xme(i,1)
                  yy = yme(i,1)
                  zz = sig(k)
               elseif (s.eq.2) then
                  xx = xme(i,2)
                  yy = yme(i,2)
                  zz = sig(k)
               elseif (s.eq.3) then
                  xx = xme(i,3)
                  yy = yme(i,3)
                  zz = sig(k)
               elseif (s.eq.4) then
                  xx = xc(i)
                  yy = yc(i)
                  if (k.eq.NZ) then
                     zz = sigv(NZ-1)+dsigv(NZ-1)
                  else
                     zz = sigv(k)
                  endif
               elseif (s.eq.5) then
                  xx = xc(i)
                  yy = yc(i)
                  if (k.eq.1) then
                     zz = sigv(1)-dsigv(1)
                  else
                     zz = sigv(k)
                  endif
               endif
!              ------------------------------------------ 
!              Approximation
               funApprox(i,k) = phi(i,k)                &
                              + dfdx(i,k)*(xx-xc(i))    &
                              + dfdy(i,k)*(yy-yc(i))    &
                              + dfdsig(i,k)*(zz-sig(k))
!              ------------------------------------------
!              Exact solution 
               if (TestDimension.eq.2) then
                  funExact(i,k) = funSolExam2D(xx,yy)
               elseif (TestDimension.eq.3) then
                  funExact(i,k) = funSolExam3D(xx,yy,zz)   
               endif
           enddo
         enddo
!        ________________________________________________________
!       |                                                        |
!       |                         Error                          |
!       |________________________________________________________|

!        --------------------------------------------------
!        Error (Maximum norm: p=infty) 
         do k=1,NZ 
            do i=1,N_CELL0 
               ErrorA1(i,k) = 0.0d0
               ErrorR1(i,k) = 0.0d0  
            enddo
         enddo

         maxErrorA1(s) = 0.0d0
         maxErrorR1(s) = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
!              -----------------
!              Absolute
               ErrorA1(i,k)= abs(funApprox(i,k)-funExact(i,k))
               maxErrorA1(s)  = max(maxErrorA1(s),ErrorA1(i,k))
!              -----------------
!              Relative
               ErrorR1(i,k) = abs(funExact(i,k)) 
               maxErrorR1(s) = max(maxErrorR1(s),ErrorR1(i,k))
            enddo
         enddo
         maxErrorR1(s) = maxErrorA1(s)/maxErrorR1(s)

!        --------------------------------------------------
!        Error (Frobenius norm: p=2)         
         sumErrorA(s) = 0.0d0
         sumErrorR(s) = 0.0d0
         do k=2,NZ-1
            do i=1,N_CELL0  
               sumErrorA(s) = sumErrorA(s) +(ErrorA1(i,k))**2
               sumErrorR(s) = sumErrorR(s) +(ErrorR1(i,k))**2
            enddo
         enddo
         sumErrorA(s) = dsqrt(sumErrorA(s)/(N_CELL0*(NZ-2)))
         sumErrorR(s) = dsqrt(sumErrorR(s)/(N_CELL0*(NZ-2)))
         sumErrorR(s) = sumErrorA(s)/sumErrorR(s) 

!        ________________________________________________________
!       |                                                        |
!       |              Solution to tecplot cell_center           |
!       |________________________________________________________|

         do k=1,NZ 
            do i=1,N_CELL0 
!              ----------------------------------------
!              Horizontal edges
               if (s.eq.1) funEdge1(i,k) = ErrorA1(i,k)
               if (s.eq.2) funEdge2(i,k) = ErrorA1(i,k)
               if (s.eq.3) funEdge3(i,k) = ErrorA1(i,k)
!              ----------------------------------------
!              Vertical edges: top (choose s=5 to bottom)
               if (s.eq.4) then
                  SolAppro(i,k) = funApprox(i,k)
                  SolExact(i,k) = funExact(i,k)
                  SolError(i,k) = ErrorA1(i,k)
               endif
            enddo
         enddo
      ENDDO    

      call SavetecEdge(funEdge1,funEdge2,funEdge3,sig)

!      ________________________________________________________
!     |                                                        |
!     |                      Display solution                  |
!     |________________________________________________________|

9     format(t10,i3,t15,e10.3,t27,e10.3,t39,e10.3,t51,e10.3)
8     format(t10,60a)
      write(*,8)'===================================================='
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'                TEST: Values at the edge            '
      write(*,'(t22,a18,i3)')'             NN = ',NN  
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'Face  Max norm    L-2 norm    Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      write(*,8)'                                                    '
      do s=1,Nface
         write(*,9) s,maxErrorA1(s),sumErrorA(s),&
                      maxErrorR1(s),sumErrorR(s)
      enddo
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

      deallocate(phi,Averagephi)
      deallocate(dfdx,dfdy,dfdsig)
      deallocate(ErrorA1,ErrorR1)
      deallocate(Exactdfdx,Exactdfdy,Exactdfdsig)
      deallocate(funEdge1,funEdge2,funEdge3)

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TEST'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
