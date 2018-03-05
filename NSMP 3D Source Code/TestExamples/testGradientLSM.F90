!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                               TESTS                                 !
!                              Feb 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testGradientLSM(SolAppro,SolExact,SolError, &
                                 xc,yc,sig,dsig,No_cp,nbe,sigv)    
   
!---------------------------------------------------------------------!
!                                                                     !
!    This program test the different subroutines of the program.      !  
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name      |    Size     | Description                     |  !  
!  |_______________|_____________|_________________________________|  !
!  | <--  fun      | (N_CELL0,NZ)| Solution of the test            |  !   
!  |_______________|_____________|_________________________________|  !
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
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      real*8,dimension(:)   :: sigv(NZ-1)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: phi
!     -----------------------------------
      real*8,dimension(:,:),allocatable  :: dfdx
      real*8,dimension(:,:),allocatable  :: dfdy
      real*8,dimension(:,:),allocatable  :: dfdsig
      real*8,dimension(:,:),allocatable  :: Exactdfdx
      real*8,dimension(:,:),allocatable  :: Exactdfdy
      real*8,dimension(:,:),allocatable  :: Exactdfdsig
!     -----------------------------------
      real*8,dimension(:,:),allocatable :: ErrorA1,ErrorA2,ErrorA3
      real*8,dimension(:,:),allocatable :: ErrorR1,ErrorR2,ErrorR3
      real*8 :: MaxErrorA1,MaxErrorA2,MaxErrorA3
      real*8 :: MaxErrorR1,MaxErrorR2,MaxErrorR3
      integer:: DisplayThis,elem
!     -----------------------------------
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: Neumanndfdn2D,Neumanndfdn3D
      real*8 :: x,y,z,fB,dfBdn,nnx,nny,nnz
!     -----------------------------------
      real*8 :: dfdxExam2D,dfdyExam2D
      real*8 :: dfdxExam3D,dfdyExam3D,dfdzExam3D
!     -----------------------------------
      integer, parameter :: DisplayThis2D = 1
      integer, parameter :: DisplayThis3D = 1

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,*) ''
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: TEST GradientLSM'
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


      allocate(phi(N_CELL,NZ))

      allocate(dfdx(N_CELL,NZ),&
               dfdy(N_CELL,NZ),&
               dfdsig(N_CELL,NZ))

      allocate(Exactdfdx(N_CELL,NZ),&
               Exactdfdy(N_CELL,NZ),&
               Exactdfdsig(N_CELL,NZ))

      allocate(ErrorA1(N_CELL0,NZ),ErrorR1(N_CELL0,NZ))
      allocate(ErrorA2(N_CELL0,NZ),ErrorR2(N_CELL0,NZ))
      allocate(ErrorA3(N_CELL0,NZ),ErrorR3(N_CELL0,NZ))

!*********************************************************************!
!                                                                     !
!                            2D Functions                             !
!                                                                     !
!*********************************************************************!

      IF (DisplayThis2D.eq.1) THEN

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
                         x = xe(i,j)
                         y = ye(i,j)
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
!                            3D Functions                             !
!                                                                     !
!*********************************************************************!

      IF (DisplayThis3D.eq.1) THEN
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
!                       Exact gradient solutions                      !
!                                                                     !
!*********************************************************************!

      IF (DisplayThis2D.eq.1) THEN
         do k=1,NZ 
            do i=1,N_CELL0
               x = xc(i)
               y = yc(i)
               Exactdfdx(i,k)   = dfdxExam2D(x,y)
               Exactdfdy(i,k)   = dfdyExam2D(x,y)
               Exactdfdsig(i,k) = 0.0d0
            enddo
         enddo
      ENDIF      
      IF (DisplayThis3D.eq.1)  THEN    
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
!                    Approximation of the Gradient                    !
!                                                                     !
!*********************************************************************!

!     _____________________________________________________
!     Tags nbe(i) for this particular example 

!     ----------------------------------
!     Dirichlet
      if (ChooseBoundary.eq.1) then
         write(*,'(t9,60a)')' WARNING!!!: All tags has been changed to nbe(i)=1'
         do i=1,N_CELL0
            nbe(i) = 1
         enddo
!     ----------------------------------
!     Dirichlet
      elseif (ChooseBoundary.eq.1) then
         write(*,'(t9,60a)')' WARNING!!!: All tags has been changed to nbe(i)=2'
         do i=1,N_CELL0
            nbe(i) = 2
         enddo
      endif

!     _____________________________________________________
!     nbe(i) tags for this particular example 

      call grandientLSM(dfdx,dfdy,dfdsig,phi,No_cp,nbe,sig) 
      !call grandientGGF(dfdx,dfdy,dfdsig,phi,xc,yc,sig,No_cp,nbe) 


!*********************************************************************!
!                                                                     !
!                               Error                                 !
!                                                                     !
!*********************************************************************!
    
!     _____________________________________________________
!     Calculate Error
      maxErrorA1 = 0.0d0
      maxErrorA2 = 0.0d0
      maxErrorA3 = 0.0d0
      maxErrorR1 = 0.0d0
      maxErrorR2 = 0.0d0
      maxErrorR3 = 0.0d0
      do k=1,NZ 
         do i=1,N_CELL0 
!           --------------------------------------------
!           Absolute
            ErrorA1(i,k) = abs(Exactdfdx(i,k)-dfdx(i,k))
            ErrorA2(i,k) = abs(Exactdfdy(i,k)-dfdy(i,k))
            ErrorA3(i,k) = abs(Exactdfdsig(i,k)-dfdsig(i,k))
            maxErrorA1   = max(maxErrorA1,ErrorA1(i,k))
            maxErrorA2   = max(maxErrorA2,ErrorA2(i,k))
            maxErrorA3   = max(maxErrorA3,ErrorA3(i,k))
!           --------------------------------------------
!           Relative
            ErrorR1(i,k) = abs(Exactdfdx(i,k)) 
            ErrorR2(i,k) = abs(Exactdfdy(i,k)) 
            ErrorR3(i,k) = abs(Exactdfdsig(i,k)) 
            maxErrorR1 = max(maxErrorR1,ErrorR1(i,k))
            maxErrorR2 = max(maxErrorR2,ErrorR2(i,k))
            maxErrorR3 = max(maxErrorR3,ErrorR3(i,k))
         enddo
      enddo
      maxErrorR1 = maxErrorA1/maxErrorR1
      maxErrorR2 = maxErrorA2/maxErrorR2
      maxErrorR3 = maxErrorA3/maxErrorR3

!     _____________________________________________________
!     Display error

8     format(t10,60a)
7     format(t18,a22,e10.3)
      write(*,8)'===================================================='
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'             TEST: GRADIENT cell-centers            '
      write(*,'(t30,a4,i3)') ' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,7)'Absolute Error dfdx = ',maxErrorA1
      write(*,7)'Absolute Error dfdy = ',maxErrorA2
      write(*,7)'Absolute Error dfdz = ',maxErrorA3
      write(*,*) ''
      write(*,7)'Relative Error dfdx = ',maxErrorR1
      write(*,7)'Relative Error dfdy = ',maxErrorR2
      write(*,7)'Relative Error dfdz = ',maxErrorR3
      write(*,8)'----------------------------------------------------'
      write(*,8)'===================================================='
      write(*,*) ''

!     _____________________________________________________
!     Save variable to display
      do k=1,NZ 
         do i=1,N_CELL0
            SolAppro(i,k) = dfdx(i,k)
            SolExact(i,k) = Exactdfdx(i,k)
            SolError(i,k) = ErrorA1(i,k)     
         enddo
      enddo


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(phi)
      deallocate(dfdx,dfdy,dfdsig)
      deallocate(Exactdfdx,Exactdfdy,Exactdfdsig)
      deallocate(ErrorA1,ErrorR1)
      deallocate(ErrorA2,ErrorR2)
      deallocate(ErrorA3,ErrorR3)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '<<<<< End   subroutine: TEST Gradient LSM'
      write(*,*) ' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
