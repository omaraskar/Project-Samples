!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM             !
!                             May 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AdvDiffEqn3D(phiNew,phivNew,                   &
                              uu,vv,ww,Gamx,Gamy,Gamz,rhs,      &  
                              phi,xc,yc,sig,dsig,No_cp,nbe,     &
                              phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                              tagBC)              
 
!---------------------------------------------------------------------!   
!                                                                     !
!     This subroutine calculates the 3D advection-diffusion equation  !
!     given the velocity profile: (uu,vv,ww),the diffusive components:!
!     (Gamx,Gamy,Gamy) and a right-hand side: rhs as follows:         !
!                                                                     !
!     d(phi)/dt  + ADV(uu,vv,ww,phi) = Diff(Gamx,Gamy,Gamz,phi) + rhs !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !  
!  |____________|____________|_____________________________________|  ! 
!  | <--phiNew  |(N_CELL,NZ) | Cell-center solution at t(n+1)      |  !
!  | <--phivNew |(N_VERT,NZ) | Cell-vertex solution at t(n+1)      |  ! 
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> uu     |(N_CELL,NZ) | x-velocity component                |  !
!  | --> vv     |(N_CELL,NZ) | y-velocity component                |  !
!  | --> ww     |(N_CELL,NZ) | z-velocity component                |  !
!  | --> Gamx   |(N_CELL,NZ) | Diffusive coefficient in x          |  !
!  | --> Gamy   |(N_CELL,NZ) | Diffusive coefficient in y          |  !
!  | --> Gamz   |(N_CELL,NZ) | Diffusive coefficient in z          |  !
!  | --> rhs    |(N_CELL,NZ) | right-hand side of the problem      |  !
!  |____________|____________|_____________________________________|  !
!  | --> phi    |(N_CELL,NZ) | Cell-center solution at t(n)        |  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> phiv   |(N_VERT,NZ) | Cell-vertex solution at t(n)        |  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
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
!  |    Am0     |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  |    AmT     |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  |    AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  !
!  |____________|____________|_____________________________________|  !
!  |    AmG     |(N_CELL0,NZ)| Gradient contribution  (advection)  |  !
!  |____________|____________|_____________________________________|  ! 
!  |    Bmv1T   |(N_CELL0,NZ)| matrix coeff. vertex 1 top    (diff)|  ! 
!  |    Bmv2T   |(N_CELL0,NZ)| matrix coeff. vertex 2 top    (diff)|  ! 
!  |    Bmv3T   |(N_CELL0,NZ)| matrix coeff. vertex 3 top    (diff)|  ! 
!  |    Bmv1B   |(N_CELL0,NZ)| matrix coeff. vertex 1 bottom (diff)|  ! 
!  |    Bmv2B   |(N_CELL0,NZ)| matrix coeff. vertex 2 bottom (diff)|  ! 
!  |    Bmv3B   |(N_CELL0,NZ)| matrix coeff. vertex 3 bottom (diff)|  !
!  |____________|____________|_____________________________________|  !  
!  |    bm      |(N_CELL0,NZ)| right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - advection3D                ( advection3D.F90 )            |  !
!  |   - diffusion3D                ( diffusion3D.F90 )            |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCcellcenter3D             ( BCcellcenter3D.F90 )         |  !
!  |   - BCvertex3D                 ( BCvertex3D.F90 )             |  !
!  |   - solSOR3D                   ( NewSOR3D.F90 )               |  !
!  |   - solGMRES3D                 ( NewSOR3D.F90 )               |  !
!  |   - gmres3DAdvDiff             ( gmres3DAdvDiff.F90 )         |  !
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

      real*8,dimension(:,:) :: phiNew(N_CELL,NZ)
      real*8,dimension(:,:) :: phivNew(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     ----------------------------------------
      integer:: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmG(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      real*8,dimension(:,:) :: bm(N_CELL0,NZ)  
!     --------------------------------------
      real*8 :: Vol,som
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      integer, parameter :: ChooseSolver = 1

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: AdvDiffEqn3D'
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
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
!     ----------------------------------------
!     Boundary Conditions 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                      phi,xc,yc,sig,dsig,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |                      Initialization                    |
!     |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0	
            Am0(i,k)  = 0.0d0
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
!           -----------------
            AmG(i,k)  = 0.0d0 
!           -----------------
            Bmv1T(i,k)= 0.0d0 
            Bmv2T(i,k)= 0.0d0 
            Bmv3T(i,k)= 0.0d0 
            Bmv1B(i,k)= 0.0d0 
            Bmv2B(i,k)= 0.0d0 
            Bmv3B(i,k)= 0.0d0 
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Advection contribution                 |
!     |________________________________________________________|

      call advection3D(Am0,Am1,Am2,Am3,AmT,AmB,AmG,   &
                       uu,vv,ww,                      &
                       phi,xc,yc,sig,No_cp,nbe,       &  
                       sigv,dsigv)
!      ________________________________________________________
!     |                                                        |
!     |                 Diffusion contribution                 |
!     |________________________________________________________|

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       -Gamx,-Gamy,-Gamz,                   &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)   
!      ________________________________________________________
!     |                                                        |
!     |                        Solution                        |
!     |________________________________________________________|

!     _________________________________________________________
!     FullExplicit             

#     ifdef KeyAdvFullExplicit
!        ______________________
!        Update inside values
         do k=2,NZ-1
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New rhs
               som =  Am0(i,k)*phi(i,k)       &
                    + Am1(i,k)*phi(jc1,k)     &
                    + Am2(i,k)*phi(jc2,k)     &
                    + Am3(i,k)*phi(jc3,k)     &
                    + AmT(i,k)*phi(i,k+1)     &       
                    + AmB(i,k)*phi(i,k-1)     &
                    + AmG(i,k)                &
                    + Bmv1T(i,k)*phiv(jv1,k)  &
                    + Bmv2T(i,k)*phiv(jv2,k)  &
                    + Bmv3T(i,k)*phiv(jv3,k)  &
                    + Bmv1B(i,k)*phiv(jv1,k-1)&
                    + Bmv2B(i,k)*phiv(jv2,k-1)&
                    + Bmv3B(i,k)*phiv(jv3,k-1)  
                bm(i,k) = -som/Vol + rhs(i,k)
!               ----------------
!               Update
                phiNew(i,k) = phi(i,k)+dt*bm(i,k)
            enddo
         enddo
!        ______________________
!        Boundary conditions
         call BCcellcenter3D(phiNew,xc,yc,sig,dsig,No_cp,nbe)
#     endif

!     _________________________________________________________
!     Explicit  

#     ifdef KeyAdvExplicit
!        ______________________
!        Update inside values
         do k=2,NZ-1
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New diagonal values
               Am0(i,k)=  dt/Vol*Am0(i,k) + 1.0d0 
!              ----------------
!              New rhs
               som =  Am1(i,k)*phi(jc1,k)     &
                    + Am2(i,k)*phi(jc2,k)     &
                    + Am3(i,k)*phi(jc3,k)     &
                    + AmT(i,k)*phi(i,k+1)     &       
                    + AmB(i,k)*phi(i,k-1)     &
                    + AmG(i,k)                &
                    + Bmv1T(i,k)*phiv(jv1,k)  &
                    + Bmv2T(i,k)*phiv(jv2,k)  &
                    + Bmv3T(i,k)*phiv(jv3,k)  &
                    + Bmv1B(i,k)*phiv(jv1,k-1)&
                    + Bmv2B(i,k)*phiv(jv2,k-1)&
                    + Bmv3B(i,k)*phiv(jv3,k-1)  
               bm(i,k) = -som/Vol + rhs(i,k)
!              ----------------
!              Update
               phiNew(i,k) = (phi(i,k)+dt*bm(i,k))/Am0(i,k)
            enddo
         enddo
!        ______________________
!        Boundary conditions
         call BCcellcenter3D(phiNew,xc,yc,sig,dsig,No_cp,nbe)
#     endif

!     _________________________________________________________
!     Semi-Implicit  

#     ifdef KeyAdvSemiImplicit
!        ______________________
!        New Matrix & rhs
         do k=1,NZ
            do i=1,N_CELL0
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New matrix coeff.
               Am0(i,k)=  dt/Vol*Am0(i,k)+ 1.0d0 
               Am1(i,k) = dt/Vol*Am1(i,k)
               Am2(i,k) = dt/Vol*Am2(i,k)
               Am3(i,k) = dt/Vol*Am3(i,k)
               AmT(i,k) = dt/Vol*AmT(i,k)
               AmB(i,k) = dt/Vol*AmB(i,k)
!              ----------------
!              New rhs
               som =   AmG(i,k)                &
                     + Bmv1T(i,k)*phiv(jv1,k)  &
                     + Bmv2T(i,k)*phiv(jv2,k)  &
                     + Bmv3T(i,k)*phiv(jv3,k)  &
                     + Bmv1B(i,k)*phiv(jv1,k-1)&
                     + Bmv2B(i,k)*phiv(jv2,k-1)&
                     + Bmv3B(i,k)*phiv(jv3,k-1) 
               bm(i,k) = phi(i,k) + dt*(-som/Vol + rhs(i,k))
            enddo
         enddo
!        ______________________
!        Call the linear solver
!        ---------------
!        SOR
         if (ChooseSolver.eq.1) then 
            print*,'Function changed!! look for the orignal in old versions'
            call SOR3Dcc(phi,                          &
                          Am0,Am1,Am2,Am3,AmT,AmB,bm,   & 
                          xc,yc,sig,dsig,No_cp,nbe)
!        ---------------
!        GMRES
         elseif (ChooseSolver.eq.2) then
            print*,'Function changed!! look for the orignal in old versions'
            call GMRES3Dcc(phi,                        &
                          Am0,Am1,Am2,Am3,AmT,AmB,bm,   & 
                          xc,yc,sig,dsig,No_cp,nbe)
         endif
!        ______________________
!        Final assignation
         do k=1,NZ 
            do i=1,N_CELL
               phiNew(i,k) = phi(i,k)
            enddo
         enddo
#     endif

!     _________________________________________________________
!     Implicit by GMRES 

#     ifdef KeyAdvImplicit
!        ______________________
!        New rhs
         do k=1,NZ
            do i=1,N_CELL
               bm(i,k) = phi(i,k) + dt*rhs(i,k)
            enddo
         enddo
!        ______________________
!        Call the GMRES solver
         call gmres3DAdvDiff(phi,phiv,                    &
                             bm,uu,vv,ww,Gamx,Gamy,Gamz,  & 
                             xc,yc,sig,dsig,No_cp,nbe,    &
                             xv,yv,sigv,dsigv,No_vp,nbev)
!        ______________________
!        Update the solution
         do k=1,NZ
            do i=1,N_CELL
               phiNew(i,k) = phi(i,k)
            enddo
         enddo
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                   New vertex values                    |
!     |________________________________________________________|

!     ----------------------------------------
!     Interpolation
      call interpolation3D(phivNew,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phiNew,xc,yc,sig,dsig,No_cp,nbe)
!     ----------------------------------------
!     Boundary Conditions 
      call BCVertex3D(phivNew,xv,yv,sigv,dsigv,No_vp,nbev, &
                      phiNew,xc,yc,sig,dsig,No_cp,nbe)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), ' End   subroutine AdvDiffEqn3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF Advection-Diffusion 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
