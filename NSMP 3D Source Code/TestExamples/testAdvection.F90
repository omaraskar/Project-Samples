!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Feb 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testAdvection(SolAppro,SolExact,SolError,      &
                               SolApprov,SolExactv,SolErrorv,   & 
                               xc,yc,sig,dsig,No_cp,nbe,        &
                               xv,yv,sigv,dsigv,No_vp,nbev)              
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program test the advection term of a given example phi.     !  
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:,:) :: Adv(N_CELL0,NZ)
!     ----------------------------------------
      real*8,dimension(:)   :: phi2D(N_CELL)
      real*8,dimension(:)   :: phiv2D(N_VERT)
      real*8,dimension(:)   :: Am02D(N_CELL0)
      real*8,dimension(:)   :: Am12D(N_CELL0)
      real*8,dimension(:)   :: Am22D(N_CELL0)
      real*8,dimension(:)   :: Am32D(N_CELL0)
      real*8,dimension(:)   :: AmG2D(N_CELL0)  
      real*8,dimension(:)   :: uu2D(N_CELL)
      real*8,dimension(:)   :: vv2D(N_CELL)
!     ----------------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmG(N_CELL0,NZ)
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:,:) :: funApprox(N_CELL,NZ)
      real*8,dimension(:,:) :: funExact(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
!     ----------------------------------------
      real*8 :: x,y,z
      real*8 :: funSolExam2D,funSolExam3D
      real*8 :: funAdvExam2D,funAdvExam3D
      real*8 :: Volume
!     ----------------------------------------
      integer:: jc1,jc2,jc3,s,jc,jj,jv1,jv2

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,*) ''
      write(*,'(t5,60a)'), '>>>>> Begin subroutine: TEST ADVECTION'
      write(*,*) ''
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!*********************************************************************!
!                                                                     !
!                          Approximation  2D                          !
!                                                                     !
!*********************************************************************!

      IF (TestDimension.eq.2) THEN
!         ________________________________________________________
!        |                                                        |
!        |             Function: Test at cell-center              |
!        |________________________________________________________|

         do i=1,N_CELL0
            x = xc(i)
            y = yc(i)
            phi2D(i) = funSolExam2D(x,y)
         enddo
!        ----------------------------------------------
!        Boundary Condition of the cell-center points    
         call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)
!         ________________________________________________________
!        |                                                        |
!        |                     Velocity profile                   |
!        |________________________________________________________|

         do i=1,N_CELL
            uu2D(i) = 1.0d0
            vv2D(i) = 1.0d0                 
         enddo
!         ________________________________________________________
!        |                                                        |
!        |                 Advection contribution                 |
!        |________________________________________________________|

         do i=1,N_CELL0	
            Am02D(i) = 0.0d0
            Am12D(i) = 0.0d0
            Am22D(i) = 0.0d0
            Am32D(i) = 0.0d0
            AmG2D(i) = 0.0d0 
         enddo

         call advection2D(Am02D,Am12D,Am22D,Am32D,AmG2D,   &
                          uu2D,vv2D,                       &
                          phi2D,xc,yc,No_cp,nbe)  
!         ________________________________________________________
!        |                                                        |
!        |                        Solution                        |
!        |________________________________________________________|

         do i=1,N_CELL0	
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            do k=1,NZ
	       Adv(i,k)  =  Am02D(i)*phi2D(i)     &
                          + Am12D(i)*phi2D(jc1)   &
                          + Am22D(i)*phi2D(jc2)   &
                          + Am32D(i)*phi2D(jc3)   &
                          + AmG2D(i)
            enddo
         enddo

!*********************************************************************!
!                                                                     !
!                          Approximation 3D                           !
!                                                                     !
!*********************************************************************!

      ELSEIF (TestDimension.eq.3) THEN  
!         ________________________________________________________
!        |                                                        |
!        |             Function: Test at cell-center              |
!        |________________________________________________________|

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
!         ________________________________________________________
!        |                                                        |
!        |                     Velocity profile                   |
!        |________________________________________________________|

         do k=1,NZ
            do i=1,N_CELL	
               uu(i,k) = 1.0d0
               vv(i,k) = 1.0d0
               ww(i,k) = 1.0d0                    
            enddo
         enddo
!         ________________________________________________________
!        |                                                        |
!        |                 Advection contribution                 |
!        |________________________________________________________|

         do k=1,NZ
            do i=1,N_CELL0	
               Am0(i,k) = 0.0d0
               Am1(i,k) = 0.0d0
               Am2(i,k) = 0.0d0
               Am3(i,k) = 0.0d0
               AmT(i,k) = 0.0d0
               AmB(i,k) = 0.0d0
               AmG(i,k) = 0.0d0 
           enddo
         enddo

         call advection3D(Am0,Am1,Am2,Am3,AmT,AmB,AmG,  &
                          uu,vv,ww,                     &
                          phi,xc,yc,sig,No_cp,nbe,      &  
                          sigv,dsigv)
!         ________________________________________________________
!        |                                                        |
!        |                        Solution                        |
!        |________________________________________________________|

         do k=2,NZ-1
            do i=1,N_CELL0
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
	       Adv(i,k)  =  Am0(i,k)*phi(i,k)     &
                          + Am1(i,k)*phi(jc1,k)   &
                          + Am2(i,k)*phi(jc2,k)   &
                          + Am3(i,k)*phi(jc3,k)   &
                          + AmT(i,k)*phi(i,k+1)   &       
                          + AmB(i,k)*phi(i,k-1)   &
                          + AmG(i,k)
            enddo
         enddo

      ENDIF

!*********************************************************************!
!                                                                     !
!                                 Error                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Exact solution: Advection term             |
!     |________________________________________________________|

!     _________________________________________________________
!     2D
      IF (TestDimension.eq.2) THEN
         do i=1,N_CELL
           do k=1,NZ 
              x = xc(i)
              y = yc(i)
              funExact(i,k) = funAdvExam2D(x,y)
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
                 funExact(i,k) = funAdvExam3D(x,y,z)
                 funExact(i,k) = funExact(i,k)*AreaCell(i)*dsigv(1) 
              else
                 funExact(i,k) = funAdvExam3D(x,y,z)
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
               funApprox(i,k) = Adv(i,k)
            enddo
         enddo
!     _________________________________________________________
!     3D
      ELSEIF (TestDimension.eq.3) THEN 
         do k=2,NZ-1 
            do i=1,N_CELL 
               funApprox(i,k) = Adv(i,k)
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
      do k=2,NZ-1 
         do i=1,N_CELL0 
            ErrorA(i,k) = abs(funApprox(i,k)-funExact(i,k))
            ErrorR(i,k) = abs(funExact(i,k))
            maxErrorA = max(maxErrorA,ErrorA(i,k))
            maxErrorR = max(maxErrorR,ErrorR(i,k))
            sumErrorA = sumErrorA+ErrorA(i,k)**2
            sumErrorR = sumErrorR+ErrorR(i,k)**2
         enddo
      enddo
      sumErrorA = dsqrt(sumErrorA/(N_CELL0*(NZ-2)))
      sumErrorR = dsqrt(sumErrorR/(N_CELL0*(NZ-2)))
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
      write(*,8)'               TEST: Advection term                 '
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
!                         Save tecplot solutions                      !
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
    	       phiv2D(nv) = funAdvExam2D(x,y)
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
    	          SolApprov(nv,k) = funAdvExam3D(x,y,z)  
               enddo
	    endif
!           __________
!           Vertical 
            z = sigv(1)
            SolApprov(nv,1)    = funAdvExam3D(x,y,z) 
            z = sigv(NZ-1)  
            SolApprov(nv,NZ-1) = funAdvExam3D(x,y,z) 
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
              SolExactv(nv,k) = funAdvExam2D(x,y)
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
              SolExactv(nv,k) = funAdvExam3D(x,y,z)
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
    	          !SolErrorv(nv,k) = 0.0d0
               enddo
	    endif
!           __________
!           Vertical 
            !SolErrorv(nv,1)    = 0.0d0
            !SolErrorv(nv,NZ-1) = 0.0d0
         enddo
      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(t5,60a)'), '<<<<< End   subroutine: TEST ADVECTION'
      write(*,*) ' '
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF testAdvection                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
