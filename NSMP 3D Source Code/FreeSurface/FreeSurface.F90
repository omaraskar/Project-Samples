!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     THE FREE SURFACE EQUATION                       !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FreeSurface(Hpr,h,eta,                    &
                             Hprv,hv,etav,                 &
                             uu,vv,                        & 
                             xc,yc,sig,dsig,No_cp,nbe,     &
                             xv,yv,sigv,dsigv,No_vp,nbev,  &
                             nRK)               
!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the free surface H = h + eta using the      !
!    equation:                                                        !
!                dH/dt +  d(H*qu)/dx + d(H*qv)/dy = 0                 ! 
!    where qu and qv are the averate velocity of u and v in the       !
!    sigma direction respectively.                                    !
!    We consider two options:                                         !
!    -------------------------------------------------------------    !
!    ChooseFreeMethod = 1: the function H in the advection part is    !
!    consider at the previous value                                   !
!                   dHnew/dt + ADV(Hold,qu,qv) = 0                    ! 
!                                                                     !
!    and solve everything in a explicit fashion.                      !
!    -------------------------------------------------------------    !
!    ChooseFreeMethod = 2: We solve the advection equation:           !
!                   dHnew/dt + ADV(Hnew,qu,qv) = 0                    ! 
!                                                                     !
!    We have four methods (chosen at file cppdefs.h) to solve the     !
!    problem :                                                        !
!    1) FullExplicit: The ADV term is calculated at time(n),          !
!                     thus it's sent it completly in the rhs.         !
!    2) Explicit    : Only the diagonal elements of ADV are           !
!                     calculated at time(n+1).                        !
!    3) SemiImplicit: In this case the cell-center and neighbors      !
!                     are calculated at time(n+1), the gradient       !
!                     GF is calculated at (n). We use SOR to solve    !
!                     the linear system.                              !
!    4) Implicit    : In this case the cell-center, neighbors         !
!                     and the gradient GF are calculated at           !
!                     time (n+1). The GMRES method is chosen          !
!                     to solve the system.                            !     
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <-- Hnew    |  (N_CELL)   | Approximate solution cell-center  |  !   
!  | <-- Hvnew   |  (N_VERT)   | Approximate solution vertex       |  !  
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | --> Hold   |(N_CELL)     | Old solution of the equation       |  !  
!  | --> xc,yc  |(N_CELL)     | Coordinates of the cell centers    |  !
!  | --> sig    |(NZ)         | Sigma value at the cell centers    |  !
!  | --> dsig   |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | --> No_cp  |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | --> nbe    |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | --> Hvold  |(N_VERT)     | Old solution of the equation       |  !  
!  | --> xv,yv  |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | --> sigv   |(NZ-1)       | sigma of the vertex points         |  !
!  | --> dsigv  |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !  
!  | --> No_vp  |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | --> nbev   |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!  | --> nRK    | integer     | Runge-Kutta step                   |  ! 
!  |____________|_____________|____________________________________|  !
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

      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      real*8,dimension(:)   :: uu(N_CELL,NZ)
      real*8,dimension(:)   :: vv(N_CELL,NZ)
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
!     --------------------------------------
      integer :: nRK
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:)   :: Hnew(N_CELL)
      real*8,dimension(:)   :: Hvnew(N_VERT)
      real*8,dimension(:)   :: Hold(N_CELL)
      real*8,dimension(:)   :: Hvold(N_VERT) 
!     --------------------------------------

      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: AmG(N_CELL0)
      real*8,dimension(:) ::  bm(N_CELL0)  
      real*8,dimension(:) ::  qu(N_CELL)
      real*8,dimension(:) ::  qv(N_CELL)
      real*8,dimension(:) :: rhs(N_CELL)
!     ----------------------------------------
      real*8 :: dtoVol,const,Vol,som,somu,somv
      integer:: it,jc1,jc2,jc3,jc
!     ----------------------------------------
      real*8 :: sumqux,sumquy,sumqvx,sumqvy,deter
      real*8 :: dqudx,dqvdy
      real*8 :: errorsys,residu,FS_funeta
!     ------------------------------------------
      integer,parameter :: ChooseFreeMethod = 2
!     ------------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Free Surface'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                       Components of the problem                     !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                       Initial values                   |
!     |________________________________________________________|

      do i=1,N_CELL 
         Hold(i) = eta(i) + h(i)   
      enddo
!      ________________________________________________________
!     |                                                        |
!     |        Average velocity values (integrals) & rhs       |
!     |________________________________________________________|

      do i=1,N_CELL
         somu = 0.0
         somv = 0.0
         do k=1,NZ-1
            somu = somu + uu(i,k)*dsig(k)
            somv = somv + vv(i,k)*dsig(k)
         enddo
         qu(i) = somu 
         qv(i) = somv   
         rhs(i) = 0.0d0 
      enddo

!*********************************************************************!
!                                                                     !
!                          Solve fully explicit                       !
!                                                                     !
!*********************************************************************!

      IF (ChooseFreeMethod ==1) THEN
         !_______________________________________
         ! New values H*qu
         do i=1,N_CELL
            qu(i) = Hold(i)*qu(i) 
            qv(i) = Hold(i)*qv(i)   
         enddo
         !_______________________________________
         ! Solution
         do i=1,N_CELL0 	
           !--------------------      
           ! Derivatives by LSM
	    sumqux = 0.0d0
	    sumquy = 0.0d0
	    sumqvx = 0.0d0
	    sumqvy = 0.0d0
	    do j=1,3
	       jc = No_cp(i,j)                
	       sumqux = sumqux + dxCC(i,j)*(qu(jc)-qu(i))
	       sumquy = sumquy + dyCC(i,j)*(qu(jc)-qu(i))
	       sumqvx = sumqvx + dxCC(i,j)*(qv(jc)-qv(i))
	       sumqvy = sumqvy + dyCC(i,j)*(qv(jc)-qv(i))
	    enddo
	    deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
            dqudx = (sum_yc2(i)*sumqux-sum_xcyc(i)*sumquy)/deter
	    dqvdy = (sum_xc2(i)*sumqvy-sum_xcyc(i)*sumqvx)/deter
           !-------------------      
           ! Update
            Hnew(i) = Hold(i)-dt*(dqudx+dqvdy)
         enddo
!        ______________________
!        Boundary conditions
         call FS_BCeta(Hnew,h,xc,yc,No_cp,nbe)


!*********************************************************************!
!                                                                     !
!                    Solve a 2D advection equation                    !
!                                                                     !
!*********************************************************************!

      ELSEIF (ChooseFreeMethod ==2) THEN

!         ________________________________________________________
!        |                                                        |
!        |                 Advection contribution                 |
!        |________________________________________________________|

         do i=1,N_CELL0	
            Am0(i) = 0.0d0
            Am1(i) = 0.0d0
            Am2(i) = 0.0d0
            Am3(i) = 0.0d0
            AmG(i) = 0.0d0 
         enddo

         call advection2D(Am0,Am1,Am2,Am3,AmG,&
                          qu,qv,Hold,xc,yc,No_cp,nbe) 
!         ________________________________________________________
!        |                                                        |
!        |                        Solution                        |
!        |________________________________________________________|

!        _________________________________________________________
!        FullExplicit            

#        ifdef KeyFreeFullExplicit
!           ______________________
!           Update inside values
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              -------------------
!              New rhs
               som =  Am0(i)*Hold(i)     &
                    + Am1(i)*Hold(jc1)   &
                    + Am2(i)*Hold(jc2)   &
                    + Am3(i)*Hold(jc3)   &
                    + AmG(i)
               bm(i) = -som/areaCell(i) + rhs(i)
!              -------------------
!              Update
               Hnew(i) = Hold(i)+dt*bm(i)
            enddo
!           ______________________
!           Boundary conditions
            call FS_BCeta(Hnew,h,xc,yc,No_cp,nbe)
#        endif

!        _________________________________________________________
!        Explicit 

#        ifdef KeyFreeExplicit
!           ______________________
!           Update inside values
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              -------------------
!              New diagonal values
               Am0(i)= 1.0d0 + dt/areaCell(i)*Am0(i)
!              -------------------
!              New rhs
               som =  Am1(i)*Hold(jc1)   &
                    + Am2(i)*Hold(jc2)   &
                    + Am3(i)*Hold(jc3)   &
                    + AmG(i)
               bm(i) = -som/areaCell(i) + rhs(i)
!              -------------------
!              Update
               Hnew(i) = (Hold(i)+dt*bm(i))/Am0(i)
            enddo
!           ______________________
!           Boundary conditions
            call FS_BCeta(Hnew,h,xc,yc,No_cp,nbe)
#        endif
!        _________________________________________________________
!        Semi-Implicit 

#        ifdef KeyFreeSemiImplicit
!           ______________________
!           New Matrix & rhs
            do i=1,N_CELL0 	
               Am0(i) = dt/areaCell(i)*Am0(i) + 1.0d0
               Am1(i) = dt/areaCell(i)*Am1(i)
               Am2(i) = dt/areaCell(i)*Am2(i)
               Am3(i) = dt/areaCell(i)*Am3(i)
               bm(i)  = Hold(i)+dt*(-AmG(i)/areaCell(i) + rhs(i))
            enddo
!           ______________________
!           Initial guess
            do i=1,N_CELL
               Hnew(i) = Hold(i)
            enddo
!           ______________________
!           Call the linear solver
            it=0
111         continue
            it=it+1
!           ---------------
!           Update
            errorsys = 0.0d0
            do i=1,N_CELL0
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
	       som =  Am1(i)*Hnew(jc1)&
                    + Am2(i)*Hnew(jc2)&
                    + Am3(i)*Hnew(jc3)
	       residu = (bm(i)-som)/Am0(i)-Hnew(i)
	       errorsys = errorsys + abs(residu)
	       Hnew(i) = Hnew(i) + relaxSOR*residu
            enddo
!           ---------------
!           Free surface BC 
            call FS_BCeta(Hnew,h,xc,yc,No_cp,nbe)
!           ---------------
!           Stop criteria   
            if (errorsys.lt.eps) then
               write(*,6) ' Solution S.0.R.: iters =',it,', error =',errorsys
            elseif (errorsys.gt.1.0d5) then
               write(*,6) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
               stop
            elseif(it.gt.MaxIters) then
               write(*,6) ' Non-convergence: iters =',it,', error =',errorsys
            else
               goto 111
            endif
!           ---------------
119         continue

6           format(t10,a26,i5,a9,e10.3)

#        endif

!        _________________________________________________________
!        Implicit by GMRES 

#        ifdef KeyFreeImplicit
         print*,'------------------------------------------------------'
         print*,' Warning!!!!!!! Not developed yet. We need to do some '
         print*,' changes in the program gmres2DAdvDiff.F90, specially '
         print*,' in the boundary conditions. Please call the new      '
         print*,' program: FS_gmres2DAdvDiff.F90                        '
         print*,'------------------------------------------------------'
         stop
!           ______________________
!           New rhs
            do i=1,N_CELL0 	
               rhs(i)  = Hold(i)+dt*rhs(i)
            enddo
!           ______________________
!           Call the GMRES solver
            call gmres2Dadvection(Hold,rhs,qu,qv, & 
                                  xc,yc,No_cp,nbe) 
!           ______________________
!           Update the solution
            do i=1,N_CELL
               Hnew(i) = Hold(i)
            enddo
#        endif

      ENDIF

!*********************************************************************!
!                                                                     !
!                       Vertex Approximation                          !
!                                                                     !
!*********************************************************************!

!     --------------------------------
!     Interpolation
      call interpolation2D(Hvnew,xv,yv,No_vp,nbev, &
                           Hnew,xc,yc,No_cp,nbe)
!     --------------------------------
!     BC
      do nv=1,N_VERT
         if (nbev(nv).ne.0) then
            Hvnew(nv) = FS_funeta(xv(nv),yv(nv),time) + hv(nv)
         endif
      enddo

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                      Final solution                    |
!     |________________________________________________________|

!     -------------------------
!     cell center 
      do i=1,N_CELL 
         Hpr(i) = Hnew(i)
         eta(i) = Hpr(i) - h(i)  
      enddo
!     -------------------------
!     vertex
      do nv=1,N_VERT 
         Hprv(nv) = Hvnew(nv)
         etav(nv) = Hprv(nv) - hv(nv) 
      enddo


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Free surface'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF Free sureface                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
