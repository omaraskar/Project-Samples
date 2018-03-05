!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        ERROR TIME EQUATION                      !
!                              Nov 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeError(U,UE,ErrorU,&
                               Uv,UEv,ErrorUv,             &  
                               xc,yc,sig,dsig,No_cp,nbe,   &
                               xv,yv,sigv,dsigv,No_vp,nbev)    

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the exact solution and errors of the     !
!    time problems. An important remark is that we display more       !
!    varaibles than we actually have.                                 !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- ufnpE  |(N_CELL,NZ)  | u Exact solution cell-center      |  !  
!  | <--- Erruf  |(N_CELL,NZ)  | u Error cell-center               |  !    
!  | <--- ufv    |(N_VERT,NZ-1)| u Approximate solution vertex     |  !
!  | <--- ufvE   |(N_VERT,NZ-1)| u Exact solution vertex           |  ! 
!  | <--- Errufv |(N_VERT,NZ-1)| u Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> ufnp  |(N_CELL,NZ)  | u Old solution of the equation     |  !
!  |____________|_____________|____________________________________|  !  
!  | ---> xc,yc |(N_CELL)     | Coordinates of the cell centers    |  !
!  | ---> sig   |(NZ)         | Sigma value at the cell centers    |  !
!  | ---> dsig  |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | ---> No_cp |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | ---> nbe   |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | ---> sigv  |(NZ-1)       | sigma of the vertex points         |  !
!  | ---> dsigv |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !  
!  | ---> No_vp |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | ---> nbev  |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |   - interpolation2D             ( interpolation2D.F90 )       |  !  
!  |   - interpolation3D             ( interpolation3D.F90 )       |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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

      real*8,dimension(:,:) :: U(N_CELL,NZ)
      real*8,dimension(:,:) :: UE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorU(N_CELL,NZ) 
      real*8,dimension(:,:) :: Uv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: UEv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: ErrorUv(N_VERT,NZ-1)
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:,:) :: hola(N_CELL,NZ)
      real*8,dimension(:,:) :: holav(N_VERT,NZ-1)

      real*8,dimension(:)   :: Err2D(N_CELL)
      real*8,dimension(:)   :: Err2Dv(N_VERT)
!     -------------------------------------
      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phiE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8,dimension(:) :: SaveErrorMax(4)
      real*8,dimension(:) :: SaveErrorSum(4)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
      real*8 :: Peakvalue,PeakvalueEx,PeakError
!     --------------------------------------
      real*8 :: x,y,z
      real*8 :: TimeExample2D,TimeExample3D     
      real*8 :: funExamNSu,funExamNSv,funExamNSw
      real*8 :: funExamNSp,funExamNSrhsp 
      integer:: DisplayThis,s
!     --------------------------------------


!*********************************************************************!
!                                                                     !
!                             Initialization                          !
!                                                                     !
!*********************************************************************!

      DisplayThis = 1
      if (DisplayThis.eq.1) then
!        _________________________________________________________
!        Peak value
         Peakvalue = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               Peakvalue   = max(Peakvalue,abs(U(i,k)))
            enddo
         enddo
         write(*,5) 'N =',NN,':  Peak value= ',Peakvalue
         write(*,*) ' '
!        _________________________________________________________
      endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TestTimeError'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF ((time.ge.tInisave).and.&
          (1d-08.ge.dtsave-(time-LastTimeSave))) THEN

!*********************************************************************!
!                                                                     !
!               Error of RK2, AdvEqn, DiffEqn, AdvDiffEqn             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |              Exact solution & Error 2D                 |
!     |________________________________________________________|

      IF (TestDimension.eq.2) THEN
!        ______________________________________________________
!        Cell-center 
         do k=1,NZ 
            do i=1,N_CELL
               x = xc(i)
               y = yc(i)
!              -----------------------------------
!              Exact values 
               UE(i,k) = TimeExample2D(x,y,time)
!              -----------------------------------
!              Error 
               ErrorU(i,k) = abs(U(i,k)-UE(i,k))
            enddo
         enddo
!        ______________________________________________________
!        Vertex
!        -----------------------------------
!        Exact values 
         do k=1,NZ-1 
            do nv=1,N_VERT
               x = xv(nv)
               y = yv(nv)
               UEv(nv,k) = TimeExample2D(x,y,time)
            enddo
         enddo
!        -----------------------------------
!        Error (interpolation)
         do i=1,N_CELL
            Err2D(i) = ErrorU(i,3)
         enddo
         call interpolation2D(Err2Dv,xv,yv,No_vp,nbev, &
                              Err2D,xc,yc,No_cp,nbe)
         do k=1,NZ-1  
            do nv=1,N_VERT
               ErrorUv(nv,k) = Err2Dv(nv)
            enddo
         enddo

!      ________________________________________________________
!     |                                                        |
!     |              Exact solution & Error 3D                 |
!     |________________________________________________________|

      ELSEIF (TestDimension.eq.3) THEN 
!        ______________________________________________________
!        Cell-center 
         do k=1,NZ 
            do i=1,N_CELL
               x = xc(i)
               y = yc(i)
               z = sig(k)
!              -----------------------------------
!              Exact values 
               UE(i,k) = TimeExample3D(x,y,z,time)
!              -----------------------------------
!              Error 
               ErrorU(i,k) = abs(U(i,k)-UE(i,k))
            enddo
         enddo
!        ______________________________________________________
!        Vertex
!        -----------------------------------
!        Exact values 
         do k=1,NZ-1 
            do nv=1,N_VERT
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)
               UEv(nv,k) = TimeExample3D(x,y,z,time)
            enddo
         enddo
!        -----------------------------------
!        Error (interpolation)
         call interpolation3D(ErrorUv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              ErrorU,xc,yc,sig,dsig,No_cp,nbe)
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |       Calculation of norm errors & Peak values         |
!     |________________________________________________________|

      do k=2,NZ-1 
         do i=1,N_CELL0 
            phiA(i,k) = U(i,k)
            phiE(i,k) = UE(i,k)
         enddo
      enddo
!     _________________________________________________________
!     Error
      maxErrorA = 0.0d0
      maxErrorR = 0.0d0
      sumErrorA = 0.0d0
      sumErrorR = 0.0d0
      do k=2,NZ-1 
         do i=1,N_CELL0 
            ErrorA(i,k) = abs(phiA(i,k)-phiE(i,k))
            ErrorR(i,k) = abs(phiE(i,k))
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
!     _________________________________________________________
!     Peak value

      Peakvalue = 0.0d0
      PeakvalueEx = 0.0d0
      do k=2,NZ-1 
         do i=1,N_CELL0 
            Peakvalue   = max(Peakvalue  ,abs(phiA(i,k)))
            PeakvalueEx = max(PeakvalueEx,abs(phiE(i,k)))
         enddo
      enddo
      PeakError = abs(Peakvalue-PeakvalueEx)

!      ________________________________________________________
!     |                                                        |
!     |                      Display solution                  |
!     |________________________________________________________|

      DisplayThis = 1
      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,8)'                                                    '
      if (RUNtestRK2approx.eq.1) &
      write(*,8)'            TEST: ODE Equation by RK-2              '
      if (RUNtestAdvEqn.eq.1) &
      write(*,8)'              TEST: Advection Equation              '
      if (RUNtestDiffEqn.eq.1) &
      write(*,8)'              TEST: Diffusion Equation              '
      if (RUNtestAdvDiffEqn.eq.1) &
      write(*,8)'         TEST: Advection-Diffusion Equation         '

      write(*,6)' Time = ',time
      write(*,7)' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      write(*,9) '   ',maxErrorA,sumErrorA,maxErrorR,sumErrorR
      write(*,8)'____________________________________________________'
      write(*,*) '                                                   '
      write(*,'(t15,a17,f10.5)') ' Peak: Approx = ',PeakValue  
      write(*,'(t15,a17,f10.5)') '       Exact  = ',PeakValueEx
      write(*,'(t15,a17,f10.5)') '       Error  = ',PeakError
      write(*,8)'===================================================='
      write(*,*)' '
      endif
!      ________________________________________________________
!     |                                                        |
!     |                      Save Errors                       |
!     |________________________________________________________|

       write(8100,*) time,PeakValue,PeakValueEx,PeakError,&
                     maxErrorA,sumErrorA,maxErrorR,sumErrorR
                         
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      5 format(t12,a4,i3,a18,e10.3)
      6  format(t22,a7,f8.4)
      7  format(t25,a4,i3)
      8  format(t5,60a)
      9  format(t5,a3,t9,e10.3,t21,e10.3,t34,e10.3,t46,e10.3)

      ENDIF

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TestTimeError'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  ERROR OF THE NAVIER-STOKES EQUATION                !
!                              Nov 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE testTimeErrorNS(ufnp,ufnpE,Erruf,ufv,ufvE,Errufv,    &
                                 vfnp,vfnpE,Errvf,vfv,vfvE,Errvfv,    &
                                 wfnp,wfnpE,Errwf,wfv,wfvE,Errwfv,    &
                                 pfnp,pfnpE,Errpf,pfv,pfvE,Errpfv,    &  
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 xv,yv,sigv,dsigv,No_vp,nbev)               

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the exact solution and errors of the     !
!    time problems. An important remark is that we display more       !
!    varaibles than we actually have.                                 !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- ufnpE  |(N_CELL,NZ)  | u Exact solution cell-center      |  !  
!  | <--- Erruf  |(N_CELL,NZ)  | u Error cell-center               |  !    
!  | <--- ufv    |(N_VERT,NZ-1)| u Approximate solution vertex     |  !
!  | <--- ufvE   |(N_VERT,NZ-1)| u Exact solution vertex           |  ! 
!  | <--- Errufv |(N_VERT,NZ-1)| u Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- vfnpE  |(N_CELL,NZ)  | v Exact solution cell-center      |  !  
!  | <--- Errvf  |(N_CELL,NZ)  | v Error cell-center               |  !    
!  | <--- vfv    |(N_VERT,NZ-1)| v Approximate solution vertex     |  !
!  | <--- vfvE   |(N_VERT,NZ-1)| v Exact solution vertex           |  ! 
!  | <--- Errvfv |(N_VERT,NZ-1)| v Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- wfnpE  |(N_CELL,NZ)  | w Exact solution cell-center      |  !  
!  | <--- Errwf  |(N_CELL,NZ)  | w Error cell-center               |  !    
!  | <--- wfv    |(N_VERT,NZ-1)| w Approximate solution vertex     |  !
!  | <--- wfvE   |(N_VERT,NZ-1)| w Exact solution vertex           |  ! 
!  | <--- Errwfv |(N_VERT,NZ-1)| w Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- pfnpE  |(N_CELL,NZ)  | p Exact solution cell-center      |  !  
!  | <--- Errpf  |(N_CELL,NZ)  | p Error cell-center               |  !    
!  | <--- pfv    |(N_VERT,NZ-1)| p Approximate solution vertex     |  !
!  | <--- pfvE   |(N_VERT,NZ-1)| p Exact solution vertex           |  ! 
!  | <--- Errpfv |(N_VERT,NZ-1)| p Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> ufnp  |(N_CELL,NZ)  | u Old solution of the equation     |  !
!  | ---> vfnp  |(N_CELL,NZ)  | v Old solution of the equation     |  !  
!  | ---> wfnp  |(N_CELL,NZ)  | w Old solution of the equation     |  !  
!  | ---> pfnp  |(N_CELL,NZ)  | p Old solution of the equation     |  !    
!  | ---> xc,yc |(N_CELL)     | Coordinates of the cell centers    |  !
!  | ---> sig   |(NZ)         | Sigma value at the cell centers    |  !
!  | ---> dsig  |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | ---> No_cp |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | ---> nbe   |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | ---> sigv  |(NZ-1)       | sigma of the vertex points         |  !
!  | ---> dsigv |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !  
!  | ---> No_vp |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | ---> nbev  |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - interpolation3D             ( interpolation3D.F90 )       |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
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

      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: ufnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Erruf(N_CELL,NZ) 
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: ufvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errufv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errvf(N_CELL,NZ) 
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errvfv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errwf(N_CELL,NZ) 
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errwfv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errpf(N_CELL,NZ) 
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errpfv(N_VERT,NZ-1)
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:) :: Err2D(N_CELL)
      real*8,dimension(:) :: Err2Dv(N_VERT)
!     -------------------------------------
      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phiE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8,dimension(:) :: SaveErrorMax(4)
      real*8,dimension(:) :: SaveErrorSum(4)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
      real*8 :: Peakvalue,PeakvalueEx,PeakError
!     --------------------------------------
      real*8 :: x,y,z
      real*8 :: TimeExample2D,TimeExample3D     
      real*8 :: funExamNSu,funExamNSv,funExamNSw
      real*8 :: funExamNSp,funExamNSrhsp 
      integer:: DisplayThis,s
!     --------------------------------------

!*********************************************************************!
!                                                                     !
!                             Initialization                          !
!                                                                     !
!*********************************************************************!

      DisplayThis = 1
      if (DisplayThis.eq.1) then
!        _________________________________________________________
!        Peak value
         Peakvalue = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               Peakvalue   = max(Peakvalue,abs(ufnp(i,k)))
            enddo
         enddo
         write(*,5) 'N =',NN,':  Peak value= ',Peakvalue
         write(*,*) ' '
!        _________________________________________________________
      endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TestTimeError'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF ((time.ge.tInisave).and.&
          (1d-08.ge.dtsave-(time-LastTimeSave))) THEN

!*********************************************************************!
!                                                                     !
!                        Navier-Stokes solution                       !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |            Exact solution & Absolute Errors            |
!     |________________________________________________________|

!     _________________________________________________________
!     Cell-center
      do k=1,NZ 
         do i=1,N_CELL
            x = xc(i)
            y = yc(i)
            z = sig(k)
!           -----------------------------------
!           Exact values 
            ufnpE(i,k) = funExamNSu(x,y,z,time)
            vfnpE(i,k) = funExamNSv(x,y,z,time)
            wfnpE(i,k) = funExamNSw(x,y,z,time) 
            pfnpE(i,k) = funExamNSp(x,y,z,time) 
!           -----------------------------------
!           Error 
            Erruf(i,k) = abs(ufnp(i,k)-ufnpE(i,k))
            Errvf(i,k) = abs(vfnp(i,k)-vfnpE(i,k))
            Errwf(i,k) = abs(wfnp(i,k)-wfnpE(i,k))
            Errpf(i,k) = abs(pfnp(i,k)-pfnpE(i,k))
         enddo
      enddo
!     -----------------------------------
!     Error BC
      do k=1,NZ  
         do i=1,N_CELL0
            if (nbe(i).ne.0) then
               Erruf(i,k) = 0.0d0
               Errvf(i,k) = 0.0d0
               Errwf(i,k) = 0.0d0 
               Errpf(i,k) = 0.0d0
            endif
            if (k.eq.1)  Erruf(i,k) = 0.0d0
            if (k.eq.NZ) Erruf(i,k) = 0.0d0
         enddo
      enddo  

!     _________________________________________________________
!     Vertices
!     -----------------------------------
!     Approx vetex
      call interpolation3D(ufv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           ufnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(vfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           vfnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(wfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           wfnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(pfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           pfnp,xc,yc,sig,dsig,No_cp,nbe)
!     -----------------------------------
!     Exact vertex  
      do k=1,NZ-1  
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
            z = sigv(k)
            ufvE(nv,k) = funExamNSu(x,y,z,time)
            vfvE(nv,k) = funExamNSv(x,y,z,time)
            wfvE(nv,k) = funExamNSw(x,y,z,time) 
            pfvE(nv,k) = funExamNSp(x,y,z,time)
         enddo
      enddo
!     -----------------------------------
!     Error (interpolation)
      call interpolation3D(Errufv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Erruf,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Errvfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Errvf,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Errwfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Errwf,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Errpfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Errpf,xc,yc,sig,dsig,No_cp,nbe)
!     -----------------------------------
!     Error BC
      do k=1,NZ-1  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
               Errufv(nv,k) = 0.0d0
               Errvfv(nv,k) = 0.0d0
               Errwfv(nv,k) = 0.0d0 
               Errpfv(nv,k) = 0.0d0
            endif
            if (k.eq.1)    Errufv(nv,k) = 0.0d0
            if (k.eq.NZ-1) Errufv(nv,k) = 0.0d0
         enddo
      enddo      
!      ________________________________________________________
!     |                                                        |
!     |                      Norm Errors                       |
!     |________________________________________________________|

      DisplayThis = 1
      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,8)'                                                    '
      write(*,8)'           TEST: Navier-Stokes problem              '
      write(*,6)' Time = ',time
      write(*,7)' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      endif

      DO s=1,4
!        ________________________________________________________
!        Norm errors 
         maxErrorA = 0.0d0
         maxErrorR = 0.0d0
         sumErrorA = 0.0d0
         sumErrorR = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               if (s.eq.1) then
                  phiA(i,k) = ufnp(i,k)
                  phiE(i,k) = ufnpE(i,k)
               elseif (s.eq.2) then
                  phiA(i,k) = vfnp(i,k)
                  phiE(i,k) = vfnpE(i,k)
               elseif (s.eq.3) then
                  phiA(i,k) = wfnp(i,k)
                  phiE(i,k) = wfnpE(i,k)
               elseif (s.eq.4) then
                  phiA(i,k) = pfnp(i,k)
                  phiE(i,k) = pfnpE(i,k)
               endif
               ErrorA(i,k) = abs(phiA(i,k)-phiE(i,k))
               ErrorR(i,k) = abs(phiE(i,k))
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
!        ________________________________________________________
!        Save norm errors 
         SaveErrorMax(s) = maxErrorA
         SaveErrorSum(s) = sumErrorA 
!        ________________________________________________________
!        Display solution 
         if (DisplayThis.eq.1) then
            if (s.eq.1) write(*,9) ' u ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.2) write(*,9) ' v ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.3) write(*,9) ' w ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.4) write(*,9) ' p ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
         endif
      ENDDO
      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,*)' '
      endif
!      ________________________________________________________
!     |                                                        |
!     |                      Save Errors                       |
!     |________________________________________________________|

      write(8100,*) time,SaveErrorMax(1),SaveErrorSum(1),&
                         SaveErrorMax(2),SaveErrorSum(2),&
                         SaveErrorMax(3),SaveErrorSum(3),&
                         SaveErrorMax(4),SaveErrorSum(4)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      5 format(t12,a4,i3,a18,e10.3)
      6  format(t22,a7,f8.4)
      7  format(t25,a4,i3)
      8  format(t5,60a)
      9  format(t5,a3,t9,e10.3,t21,e10.3,t34,e10.3,t46,e10.3)

      ENDIF

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TestTimeError'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   END OF TEST NAVIER-STOKES EQUATION                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
