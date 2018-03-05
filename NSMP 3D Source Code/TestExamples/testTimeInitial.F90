!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INITIALIZATION OF THE MAIN VARIABLES              !
!                             Feb 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeInitial(ufn,vfn,wfn,pfn,           &
                                 ufv,vfv,wfv,pfv,           &
                                 xc,yc,sig,dsig,No_cp,nbe,  &
                                 xv,yv,sigv,dsigv,No_vp,nbev)              
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the initial condition of the time test   !  
!    problems.                                                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |           Name       |    Size     | Description              |  !  
!  |______________________|_____________|__________________________|  !
!  | <--- ufn,vfn,wfn,pfn |(N_CELL,NZ)  | Initial solution center  |  !
!  | <--- ufv,vfv,wfv,pfv |(N_VERT,NZ-1)| Initial solution vertex  |  !
!  |______________________|_____________|__________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc |(N_CELL)    | Coordinates of the cell centers     |  !
!  | ---> sig   |(NZ)        | Sigma value at the cell centers     |  !
!  | ---> dsig  |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | ---> No_cp |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | ---> nbe   |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | ---> sigv  |(NZ-1)      | sigma of the vertex points          |  !
!  | ---> dsigv |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
!  | ---> No_vp |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | ---> nbev  |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
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
!     |      Declaration of variables      |
!     |____________________________________|

      real*8,dimension(:,:) :: ufn(N_CELL,NZ)
      real*8,dimension(:,:) :: vfn(N_CELL,NZ)
      real*8,dimension(:,:) :: wfn(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
      integer,dimension(:)  :: nbev(N_VERT)  
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: AverageTestc
      real*8 :: summ,summ1,summ2
      real*8 :: x1,x2,x3,y1,y2,y3,z1,z2
      real*8 :: a11,a21,a12
      real*8 :: a22,a13,a23
      real*8 :: xmy1,xmy2,xmy3
      real*8 :: xpy1,xpy2,xpy3
      real*8 :: cp,b1,b2
      real*8 :: xx,yy,zz
      real*8 :: dtoVol
      integer:: DisplayThis
      integer:: jv1,jv2,jv3,s,jc
      integer:: jc1,jc2,jc3
      real*8 :: mx1,mx2,mx3
      real*8 :: my1,my2,my3
      real*8 :: ccx,ccy,ccz,rr,r1,r2,r3,ra
      real*8 :: fB
!     ----------------------------------------
      real*8 :: x,y,z,TimeExample2D,TimeExample3D
      real*8 :: funExamNSu,funExamNSv,funExamNSw,funExamNSp
!     ----------------------------------------
      real*8 :: Peakvalue,PeakvalueEx,PeakError
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TestTimeInitial'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Initial time of the Example                     !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Navier-Stokes Equation at t=0.0            |
!     |________________________________________________________|
 
      IF (RUNTestNSEqn.eq.1) THEN
!        ______________________________________________________
!        Initial Cell center
         do i=1,N_CELL
            do k=1,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)
               ufn(i,k) = funExamNSu(x,y,z,time)
               vfn(i,k) = funExamNSv(x,y,z,time)
               wfn(i,k) = funExamNSw(x,y,z,time)
               pfn(i,k) = funExamNSp(x,y,z,time)
            enddo
         enddo
!        ______________________________________________________
!        InitialVertex
         do nv=1,N_VERT
            do k=1,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)                                                   
    	       ufv(nv,k) = funExamNSu(x,y,z,time) 
               vfv(nv,k) = funExamNSv(x,y,z,time) 
               wfv(nv,k) = funExamNSw(x,y,z,time) 
               pfv(nv,k) = funExamNSp(x,y,z,time)
            enddo
         enddo
!        ______________________________________________________
!        Write initial time error file 

         maxErrorA = 0.0d0
         sumErrorA = 0.0d0
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
            write(8100,*) time,maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
         IF (rang_topo.eq.0) THEN
            write(8100,*) time,maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA,&
                               maxErrorA,sumErrorA
         ENDIF
#        endif
!        =============== END ================    
!        ====================================

      ELSE
!      ________________________________________________________
!     |                                                        |
!     |       Other initial examples: TimeExample at t=0       |
!     |________________________________________________________|

!        ______________________________________________________
!        Initial values cell-centers        
!        2D ------------------------------------
         if (TestDimension.eq.2) then
            do i=1,N_CELL
               do k=1,NZ 
                  x = xc(i)
                  y = yc(i)
                  ufn(i,k) = TimeExample2D(x,y,time)
               enddo
            enddo
!        3D ------------------------------------
         elseif (TestDimension.eq.3) then 
            do i=1,N_CELL
               do k=1,NZ 
                  x = xc(i)
                  y = yc(i)
                  z = sig(k)
                  ufn(i,k)  = TimeExample3D(x,y,z,time)
               enddo
            enddo
         endif
!        ______________________________________________________
!        Initial values vertex
!        2D ------------------------------------
         if (TestDimension.eq.2) then
            do nv=1,N_VERT
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  ufv(nv,k) = TimeExample2D(x,y,time)
               enddo
            enddo
!        3D ------------------------------------
         elseif (TestDimension.eq.3) then 
            do nv=1,N_VERT
               do k=1,NZ-1 
                  x = xv(nv)
                  y = yv(nv)
                  z = sig(k)
                  ufv(nv,k) = TimeExample3D(x,y,time)
               enddo
            enddo
         endif
!        ______________________________________________________
!        Initial Errors & peak values
!        ---------------------------------------
!        Initial peak value

         Peakvalue = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               Peakvalue = max(Peakvalue,abs(ufn(i,k)))
            enddo
         enddo
         PeakvalueEx = Peakvalue
         PeakError = 0.0d0
!        ---------------------------------------
!        Initial errors 

         maxErrorA = 0.0d0
         sumErrorA = 0.0d0
         maxErrorR = 0.0d0
         sumErrorR = 0.0d0

!        ---------------------------------------
!        Write initial time error file

!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
            write(8100,*) time,PeakValue,PeakValueEx,PeakError,&
                          maxErrorA,sumErrorA,maxErrorR,sumErrorR
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
         IF (rang_topo.eq.0) THEN
            write(8100,*) time,PeakValue,PeakValueEx,PeakError,&
                          maxErrorA,sumErrorA,maxErrorR,sumErrorR
         ENDIF
#        endif
!        =============== END ================ 
!        ====================================

!        ---------------------------------------
!        Display initial solution

         5 format(t12,a18,e10.3)
         6 format(t22,a7,f8.4)
         7 format(t25,a4,i3)
         8 format(t5,60a)
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
            print*,'  '
            write(*,8)'____________________________________________________'
            write(*,8)'                                                    '
            write(*,8)'               TEST: Time Equation                  '
            print*,'  '
            write(*,6)' Time = ',time
            write(*,7)' N = ',NN
            write(*,5)' Peak value = ',Peakvalue
            write(*,8)'____________________________________________________'
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
         IF (rang_topo.eq.0) THEN
            print*,'  '
            write(*,8)'===================================================='
            write(*,8)'                                                    '
            write(*,8)'          TEST: Time Equation (MPI Parallel)        '
            print*,'  '
            write(*,6)' Time = ',time
            write(*,7)' N = ',NN
            write(*,5)' Peak value = ',Peakvalue
            write(*,8)'                                                    '
            write(*,8)'===================================================='
         ENDIF
#        endif
!        =============== END ================ 
!        ====================================
      ENDIF

!*********************************************************************!
!                                                                     !
!                             Extra Funtions                          !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(AverageTestc(N_CELL,NZ))

!      ________________________________________________________
!     |                                                        |
!     |                AverageTestc = int(fun)/m(V)            |
!     |________________________________________________________|
     
      DO i=1,N_CELL0
         do k=1,NZ
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            x1 = xv(jv1)
            y1 = yv(jv1) 
            x2 = xv(jv2)
            y2 = yv(jv2)
            x3 = xv(jv3)
            y3 = yv(jv3)       
            mx1 = 0.5d0*(x1+x2)
            my1 = 0.5d0*(y1+y2)
            mx2 = 0.5d0*(x2+x3)
            my2 = 0.5d0*(y2+y3)
            mx3 = 0.5d0*(x3+x1)
            my3 = 0.5d0*(y3+y1)
            z1  = 0.5d0*(sigv(k)+sigv(k+1)) 
!           --------------
!           2D
            if (TestDimension.eq.2) then
              ccx = -0.5d0
              ccy =  0.0d0
              rr = 0.25d0
              r1 = dsqrt((mx1-ccx)**2+(my1-ccy)**2)
              r2 = dsqrt((mx2-ccx)**2+(my2-ccy)**2)
              r3 = dsqrt((mx3-ccx)**2+(my3-ccy)**2)
!           --------------
!           3D       
            elseif (TestDimension.eq.3) then
              ccx = -0.5d0*(1.0d0/dsqrt(2.0d0))
              ccy =  0.5d0*(1.0d0/dsqrt(2.0d0))
              ccz = -0.5d0*(1.0d0/dsqrt(2.0d0))
              rr = 0.25d0
              r1 = dsqrt((mx1-ccx)**2+(my1-ccy)**2+(z1-ccz)**2)
              r2 = dsqrt((mx2-ccx)**2+(my2-ccy)**2+(z1-ccz)**2)
              r3 = dsqrt((mx3-ccx)**2+(my3-ccy)**2+(z1-ccz)**2)
            endif
!           ------------------------------------------------------
!           Average integral
            summ = 0.0d0
            if (r1.le.rr) summ = summ+(1.0d0/3.0d0)*dcos(2.0d0*pi*r1)&
                                                   *dcos(2.0d0*pi*r1)
            if (r2.le.rr) summ = summ+(1.0d0/3.0d0)*dcos(2.0d0*pi*r2)&
                                                   *dcos(2.0d0*pi*r2)
            if (r3.le.rr) summ = summ+(1.0d0/3.0d0)*dcos(2.0d0*pi*r3)&
                                                   *dcos(2.0d0*pi*r3)
            AverageTestc(i,k) = summ
         enddo
      ENDDO
  
      do i=1,N_CELL0
	if (nbe(i).ne.0) then	
	   do j=1,3
              fB = 0.0d0
	      nc=No_cp(i,j)
              if (nc.gt.N_CELL0) then                 
                 do k=1,NZ                      
                    AverageTestc(nc,k) = 2.0d0*fB - AverageTestc(i,k)                   
                 enddo
              endif 
	   enddo
         endif
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

      deallocate(AverageTestc)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TestTimeInitial'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                     END OF testAdvection Initial                    !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
