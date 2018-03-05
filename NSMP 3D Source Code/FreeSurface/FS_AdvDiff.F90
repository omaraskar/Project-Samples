!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!     SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM WITH FREE SURFACE   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_AdvDiff(phiNew,phivNew,                   &
                            rhs,Gamx,Gamy,Gamz,uu,vv,ww,pfn,  &  
                            phi,xc,yc,sig,dsig,No_cp,nbe,     &
                            phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,eta,                        &
                            Hprv,hv,etav,                     &
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
!  |   - BCvelcenter3D              ( BCvelocity.F90 )             |  !
!  |   - BCvelvertex3D              ( BCvelocity.F90 )             |  !
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
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
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
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
!     ----------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
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
      real*8 :: Vol,som,errorsys,residu
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
      integer :: it,DoThis
!     --------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: FS_AdvDiff'
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
      call FS_BCvelv(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                     phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC)
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
!     |                 Diffusion contribution                 |
!     |________________________________________________________|

      DoThis = 1
      if (DoThis==1) then
      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       -Gamx,-Gamy,-Gamz,                   &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)   


      do i=1,N_CELL0	
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
         do k=1,NZ
            Am0(i,k) = Am0(i,k)/Hpr(i)  
            Am1(i,k) = Am1(i,k)/Hpr(jc1)  
            Am2(i,k) = Am2(i,k)/Hpr(jc2)  
            Am3(i,k) = Am3(i,k)/Hpr(jc3)  
            AmT(i,k) = AmT(i,k)/Hpr(i)          
            AmB(i,k) = AmB(i,k)/Hpr(i)     
            Bmv1T(i,k) = Bmv1T(i,k)/Hprv(jv1)
            Bmv2T(i,k) = Bmv2T(i,k)/Hprv(jv2)
            Bmv3T(i,k) = Bmv3T(i,k)/Hprv(jv3)
            Bmv1B(i,k) = Bmv1B(i,k)/Hprv(jv1)
            Bmv2B(i,k) = Bmv2B(i,k)/Hprv(jv2)
            Bmv3B(i,k) = Bmv3B(i,k)/Hprv(jv3)
         enddo
      enddo

      endif
!      ________________________________________________________
!     |                                                        |
!     |                 Advection contribution                 |
!     |________________________________________________________|

      !call advectionVelocity(Am0,Am1,Am2,Am3,AmT,AmB,AmG,   &
      !                       uu,vv,ww,pfn,                  &
      !                       phi,xc,yc,sig,No_cp,nbe,       &  
      !                       sigv,dsigv)

      call FS_advection(Am0,Am1,Am2,Am3,AmT,AmB,AmG,   &
                        uu,vv,ww,pfn,                  &
                        phi,xc,yc,sig,No_cp,nbe,       &  
                        sigv,dsigv)
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
         call FS_BCvelc(phiNew,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC)       
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
         call FS_BCvelc(phiNew,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC) 
#     endif

!     _________________________________________________________
!     Semi-Implicit  

#     ifdef KeyAdvSemiImplicit
!        ______________________
!        New Matrix & rhs
!        ______________________

         do k=2,NZ-1
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
!        ______________________

!        ----------------------
!        Initial guess
         do k=1,NZ
            do i=1,N_CELL  
               phiNew(i,k) = phi(i,k)
            enddo
         enddo 
!        ----------------------
!        Iterations
         it=0
111      continue
         it=it+1 
         errorsys = 0.0d0
         do k=2,NZ-1
            do i=1,N_CELL0
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
	       som =  Am1(i,k)*phiNew(jc1,k) &
                    + Am2(i,k)*phiNew(jc2,k) &
                    + Am3(i,k)*phiNew(jc3,k) &
                    + AmT(i,k)*phiNew(i,k+1) &       
                    + AmB(i,k)*phiNew(i,k-1) 
	       residu = (bm(i,k)-som)/Am0(i,k)-phiNew(i,k)
	       errorsys = errorsys + abs(residu)
	       phiNew(i,k) = phiNew(i,k) + relaxSOR*residu
            enddo
         enddo
!        ----------------------
!        Boundary conditions
         call FS_BCvelc(phiNew,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC) 
!        ----------------------
!        Convergence criteria
         if (errorsys.lt.eps) then
            !write(*,7) 'Solution S0R 3D u: iters =',it,', error =',errorsys
         elseif (errorsys.gt.1.0d5) then
            write(*,7) 'DIVERGENCE !!!! u : iters =',it,', error =',errorsys
            stop
         elseif(it.gt.MaxIters) then
            write(*,7) 'Non-convergence u: iters =',it,', error =',errorsys
         else
            goto 111
         endif
119      continue
         7 format(t10,a24,i5,a9,e10.3)
#     endif

!     _________________________________________________________
!     Implicit by GMRES 

#     ifdef KeyAdvImplicit
         print*,'------------------------------------------------------'
         print*,' Warning!!!!!!! Not developed yet. We need to do some '
         print*,' changes in the program gmres3DAdvDiff.F90, specially '
         print*,' in the boundary conditions. Please call the new      '
         print*,' program: FS_gmres3DAdvDiff.F90                        '
         print*,'------------------------------------------------------'
         stop
!        ______________________
!        New rhs
         do k=1,NZ
            do i=1,N_CELL0
               bm(i) = phi(i) + dt*rhs(i)
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
               phiNew(i) = phi(i)
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
      call FS_BCvelv(phivNew,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                     phiNew,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,eta,tagBC)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), ' End   subroutine FS_AdvDiff'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF Advection-Diffusion 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
