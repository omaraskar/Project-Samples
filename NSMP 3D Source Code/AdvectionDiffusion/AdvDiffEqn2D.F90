!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM             !
!                             May 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AdvDiffEqn2D(phiNew,phivNew,          &
                              uu,vv,Gamx,Gamy,rhs,     & 
                              phi,xc,yc,No_cp,nbe,     &
                              phiv,xv,yv,No_vp,nbev,   &
                              tagBC)              
 
!---------------------------------------------------------------------!   
!                                                                     !
!     This subroutine calculates the 2D advection-diffusion equation  !
!     given the velocity profile: (uu,vv), the diffusive components:  !
!     (Gamx,Gamy) and right-hand side: rhs as follows                 !
!                                                                     !
!        d(phi)/dt  + ADV(uu,vv,phi) = Diff(Gamx,Gamy,phi) + rhs      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !  
!  |____________|____________|_____________________________________|  ! 
!  | <--phiNew  |(N_CELL)    | Cell-center solution at t(n+1)      |  !
!  | <--phivNew |(N_VERT)    | Cell-vertex solution at t(n+1)      |  ! 
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> uu     |(N_CELL)    | x-velocity component                |  !
!  | --> vv     |(N_CELL)    | y-velocity component                |  !
!  | --> Gamx   |(N_CELL)    | Diffusive coefficient in x          |  !
!  | --> Gamy   |(N_CELL)    | Diffusive coefficient in y          |  !
!  | --> rhs    |(N_CELL)    | right-hand side of the problem      |  !
!  |____________|____________|_____________________________________|  !
!  | --> phi    |(N_CELL)    | Cell-center solution at t(n)        |  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> phiv   |(N_VERT)    | Cell-vertex solution at t(n)        |  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  ! 
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbe    |(N_CELL0)   | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!  | --> tagBC  | integer    | Tag = 1:vel.u, =2:vel.v, =3:vel.w   |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     |(N_CELL0)   | matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0)   | matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0)   | matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0)   | matrix coeff. horizontal neigborn 3 |  !
!  |____________|____________|_____________________________________|  !
!  |    AmG     |(N_CELL0)   | Gradient contribution  (advection)  |  !
!  |____________|____________|_____________________________________|  !
!  |    Bmv1    |(N_CELL0)   | matrix coeff. vertex 1 (diffusion)  |  ! 
!  |    Bmv2    |(N_CELL0)   | matrix coeff. vertex 2 (diffusion)  |  ! 
!  |    Bmv3    |(N_CELL0)   | matrix coeff. vertex 3 (diffusion)  |  !
!  |____________|____________|_____________________________________|  !  
!  |    bm      |(N_CELL0)   | right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - advection2D                ( advection2D.F90 )            |  !
!  |   - diffusion2D                ( diffusion2D.F90 )            |  !
!  |   - interpolation2D            ( interpolation2D.F90 )        |  !
!  |   - BCcellcenter2D             ( BCcellcenter2D.F90 )         |  !
!  |   - BCvertex2D                 ( BCvertex2D.F90 )             |  !
!  |   - solSOR2D                   ( NewSOR2D.F90 )               |  !
!  |   - solGMRES2D                 ( NewSOR2D.F90 )               |  !
!  |   - gmres2DAdvDiff             ( gmres2DAdvDiff.F90 )         |  !
!  |_______________________________________________________________|  !
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

      real*8,dimension(:,:)  :: phiNew(N_CELL)
      real*8,dimension(:,:)  :: phivNew(N_VERT)
!     --------------------------------------
      real*8,dimension(:)    ::  uu(N_CELL)
      real*8,dimension(:)    ::  vv(N_CELL)
      real*8,dimension(:)    :: rhs(N_CELL)
      real*8,dimension(:)    :: Gamx(N_CELL)
      real*8,dimension(:)    :: Gamy(N_CELL)
!     --------------------------------------
      real*8, dimension(:,:) :: phi(N_CELL)
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     --------------------------------------
      real*8, dimension(:,:) :: phiv(N_VERT)
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT) 
!     --------------------------------------
      integer:: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: AmG(N_CELL0)
      real*8,dimension(:) :: Bmv1(N_CELL0) 
      real*8,dimension(:) :: Bmv2(N_CELL0) 
      real*8,dimension(:) :: Bmv3(N_CELL0) 
      real*8,dimension(:) ::  bm(N_CELL0)  
!     ----------------------------------------
      real*8 :: Vol,som
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      integer,parameter :: ChooseSolver = 2
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: AdvDiffEqn2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                       Advection-Diffusion  3D                       !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                      Vertex values                     |
!     |________________________________________________________|

!     ----------------------------------------
!     Interpolation
      call interpolation2D(phiv,xv,yv,No_vp,nbev, &
                           phi,xc,yc,No_cp,nbe)
!     ----------------------------------------
!     Boundary Conditions 
      call BCVertex2D(phiv,xv,yv,No_vp,nbev, &
                      phi,xc,yc,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |                      Initialization                    |
!     |________________________________________________________|

      do i=1,N_CELL0	
         Am0(i)  = 0.0d0
         Am1(i)  = 0.0d0
         Am2(i)  = 0.0d0
         Am3(i)  = 0.0d0
!        ---------------
         AmG(i)  = 0.0d0  
!        ---------------
         Bmv1(i) = 0.0d0 
         Bmv2(i) = 0.0d0 
         Bmv3(i) = 0.0d0 
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Advection contribution                 |
!     |________________________________________________________|

      call advection2D(Am0,Am1,Am2,Am3,AmG,   &
                       uu,vv,                 &
                       phi,xc,yc,No_cp,nbe)  
!      ________________________________________________________
!     |                                                        |
!     |                 Diffusion contribution                 |
!     |________________________________________________________|

      call diffusion2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                       -Gamx,-Gamy,                    &
                       xc,yc,No_cp,nbe,                &
                       xv,yv,No_vp,nbev)
!      ________________________________________________________
!     |                                                        |
!     |                        Solution                        |
!     |________________________________________________________|

!     _________________________________________________________
!     FullExplicit            

#     ifdef KeyAdvFullExplicit
!        ______________________
!        Update inside values
         do i=1,N_CELL0 	
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Vol = areaCell(i)
!           -------------------
!           New rhs
            som =  Am0(i)*phi(i)     &
                 + Am1(i)*phi(jc1)   &
                 + Am2(i)*phi(jc2)   &
                 + Am3(i)*phi(jc3)   &
                 + AmG(i)            &
                 + Bmv1(i)*phiv(jv1) &
                 + Bmv2(i)*phiv(jv2) &
                 + Bmv3(i)*phiv(jv3) 
            bm(i) = -som/Vol + rhs(i)
!           -------------------
!           Update
            phiNew(i) = phi(i) + dt*bm(i)
         enddo
!        ______________________
!        Boundary conditions
         call BCcellcenter2D(phiNew,xc,yc,No_cp,nbe)
#     endif

!     _________________________________________________________
!     Explicit 

#     ifdef KeyAdvExplicit
!        ______________________
!        Update inside values
         do i=1,N_CELL0 	
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Vol = areaCell(i)
!           -------------------
!           New diagonal values
            Am0(i)= dt/Vol*Am0(i) + 1.0d0
!           -------------------
!           New rhs
            som =  Am1(i)*phi(jc1)   &
                 + Am2(i)*phi(jc2)   &
                 + Am3(i)*phi(jc3)   &
                 + AmG(i)            &
                 + Bmv1(i)*phiv(jv1) &
                 + Bmv2(i)*phiv(jv2) &
                 + Bmv3(i)*phiv(jv3) 
            bm(i) = -som/Vol + rhs(i)
!           -------------------
!           Update
            phiNew(i) = (phi(i)+dt*bm(i))/Am0(i)
         enddo
!        ______________________
!        Boundary conditions
         call BCcellcenter2D(phiNew,xc,yc,No_cp,nbe)
#     endif

!     _________________________________________________________
!     Semi-Implicit 

#     ifdef KeyAdvSemiImplicit
!        ______________________
!        New Matrix & rhs
         do i=1,N_CELL0 
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Vol = areaCell(i)
!           ----------------
!           New matrix coeff.
            Am0(i) = dt/Vol*Am0(i) + 1.0d0 
            Am1(i) = dt/Vol*Am1(i)
            Am2(i) = dt/Vol*Am2(i)
            Am3(i) = dt/Vol*Am3(i)
!           ----------------
!           New rhs
	    som =   AmG(i)            & 
                  + Bmv1(i)*phiv(jv1) &
                  + Bmv2(i)*phiv(jv2) &
                  + Bmv3(i)*phiv(jv3)
            bm(i) = phi(i) + dt*(-som/Vol+rhs(i))
         enddo
!        ______________________
!        Call the linear solver
!        -------------------
!        SOR
         if (ChooseSolver.eq.1) then 
            print*,'Function changed!! look for the orignal in old versions'
            ! call solSOR2D(phi,                 &
            !               Am0,Am1,Am2,Am3,bm,  & 
            !               xc,yc,No_cp,nbe)
!        -------------------
!        GMRES
         elseif (ChooseSolver.eq.2) then
            print*,'Function changed!! look for the orignal in old versions'
             !call solGMRES2D(phi,               &
             !              Am0,Am1,Am2,Am3,bm,  & 
             !              xc,yc,No_cp,nbe)
         endif
!        ______________________
!        Update the solution
         do i=1,N_CELL
            phiNew(i) = phi(i)
         enddo
#     endif

!     _________________________________________________________
!     Implicit by GMRES 

#     ifdef KeyAdvImplicit
!        ______________________
!        New rhs
         do i=1,N_CELL0
            bm(i) = phi(i) + dt*rhs(i)
         enddo
!        ______________________
!        Call the GMRES solver
         call gmres2DAdvDiff(phi,phiv,           &
                             bm,uu,vv,Gamx,Gamy, &
                             xc,yc,No_cp,nbe,    &
                             xv,yv,No_vp,nbev) 
!        ______________________
!        Update the solution
         do i=1,N_CELL
            phiNew(i) = phi(i)
         enddo
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                   New vertex values                    |
!     |________________________________________________________|

!     ----------------------------------------
!     Interpolation
      call interpolation2D(phivNew,xv,yv,No_vp,nbev, &
                           phiNew,xc,yc,No_cp,nbe)
!     ----------------------------------------
!     Boundary Conditions 
      call BCVertex2D(phivNew,xv,yv,No_vp,nbev, &
                      phiNew,xc,yc,No_cp,nbe)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'End subroutine: AdvDiff2DEqn'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!!
!                        END OF Advection-Diffusion                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
