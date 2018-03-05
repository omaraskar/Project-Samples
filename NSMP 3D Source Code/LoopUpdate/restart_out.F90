!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      WRITE RESULTS TO RESTART                       !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp, & 
                             alphasnp,usnp,vsnp,wsnp,psnp, &
                             zct)

!---------------------------------------------------------------------!
!                                                                     !
!    This program writes all the variables needed to restart a        !
!    new simulation using the final values of the previous si-        !
!    mulation.                                                        ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !  
!  |______________|___________|____________________________________|  !
!  | --> alphafnp |(N_CELL,NZ)| Fluid control volume at t(n+1)     |  !
!  | --> ufnp     |(N_CELL,NZ)| Velocity component u_f at t(n+1)   |  !
!  | --> vfnp     |(N_CELL,NZ)| Velocity component v_f at t(n+1)   |  !      
!  | --> wfnp     |(N_CELL,NZ)| Velocity component w_f at t(n+1)   |  !
!  | --> pfnp     |(N_CELL,NZ)| Pressure of the fluid at t(n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> alphasnp |(N_CELL,NZ)| Solid control volume at t(n+1)     |  !
!  | --> usnp     |(N_CELL,NZ)| Velocity component u_s at t(n+1)   |  !
!  | --> vsnp     |(N_CELL,NZ)| Velocity component v_s at t(n+1)   |  !     
!  | --> wsnp     |(N_CELL,NZ)| Velocity component w_s at t(n+1)   |  !
!  | --> psnp     |(N_CELL,NZ)| Pressure of the solid at t(n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> zct      |(N_CELL,NZ)| zc at t(n+1)                       |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  |--- N_CELL   | Total number of the cells                       |  !
!  |--- NZ       | Points in the sigma direction                   |  !
!  |    time     | Final time calculated                           |  !
!  |    dt       | Time step                                       |  !
!  | SaveCounter | Counting integer of saving files                |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  | * i,k       |  Loop counters                                  |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  |   ifile     |  Writing format                                 |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |   Keys, subroutines and common parameters              |
!     |________________________________________________________|

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

!      ________________________________________________________
!     |                                                        |
!     |           Definition of variables                      |
!     |________________________________________________________|

      real*8,dimension(:,:):: alphafnp
      real*8,dimension(:,:):: ufnp
      real*8,dimension(:,:):: vfnp
      real*8,dimension(:,:):: wfnp
      real*8,dimension(:,:):: pfnp   
      real*8,dimension(:,:):: alphasnp
      real*8,dimension(:,:):: usnp
      real*8,dimension(:,:):: vsnp
      real*8,dimension(:,:):: wsnp
      real*8,dimension(:,:):: psnp
      real*8,dimension(:,:):: zct

!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|

       integer:: ifile

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: restart_out '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|

77    format('         time = ',e12.5)
78    format('           dt = ',e12.5)
79    format('  SaveCounter = ',i8)
35    format(6e12.5)
36    format(i8)
5     format(t2,80a)

!*********************************************************************!
!                                                                     !
!                      Write file (sequential)                        !
!                                                                     !
!*********************************************************************!

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
!      ________________________________________________________
!     |                                                        |
!     |                  Open file ifile                       |
!     |________________________________________________________|

      ifile=16
      write(filerepout(11:15),'(i4.4)') SaveCounter
      open(ifile,file=filerepout,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________| 

      write(ifile,77) time
      write(ifile,78) dt
      write(ifile,79) SaveCounter
!      __________________________________
!     |                                  |
!     |         Fluid variables          |
!     |__________________________________| 

      write(ifile,5) 'alphaf'
      write(ifile,35)((alphafnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'uf'
      write(ifile,35)((ufnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'vf'
      write(ifile,35)((vfnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'wf'
      write(ifile,35)((wfnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'pf'
      write(ifile,35)((pfnp(i,k),i=1,N_CELL),k=1,NZ)
!      __________________________________
!     |                                  |
!     |         Solid variables          |
!     |__________________________________| 

      write(ifile,5) 'alphas'
      write(ifile,35)((alphasnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'us'
      write(ifile,35)((usnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'vs'
      write(ifile,35)((vsnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'ws'
      write(ifile,35)((wsnp(i,k),i=1,N_CELL),k=1,NZ)
      write(ifile,5) 'ps'
      write(ifile,35)((psnp(i,k),i=1,N_CELL),k=1,NZ)

!      __________________________________
!     |                                  |
!     |           Vertical space         |
!     |__________________________________| 

      write(ifile,5) 'zc'
      write(ifile,35)((zct(i,k),i=1,N_CELL),k=1,NZ)

!      ________________________________________________________
!     |                                                        |
!     |                      Close file                        |
!     |________________________________________________________|

      close(ifile)

#     endif
!     =============== END ================    
!     ====================================


!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: restart_out'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END FINALIZATION                            !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
