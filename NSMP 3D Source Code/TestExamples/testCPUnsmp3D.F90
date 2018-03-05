!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           ---  TEST TO COMPARE WITH MPI COMMUNICATION 3D ---        !
!                      Miguel Angel Uh Zapata                         !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

       SUBROUTINE TestCPUnsmp3D(xc,yc,nbe,No_cp)

!---------------------------------------------------------------------!
!                                                                     !
!     Programe to test definitions and comunications of the MPI       !
!     formulation. I define a function which sum its three neighbors  !
!     values and itself. I calculate the sum of all independent       !
!     elements and use an MPI command to collect the values in all    !
!     processors. The final result is compare with the global solu-   !
!     tion after a determinated number of iterations.                 !
!                                                                     !
!---------------------------------------------------------------------!  

!*********************************************************************!
!                                                                     !
!                        Definition of variables                      !
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

      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer,dimension(:,:):: No_cp(N_CELL,3)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8 ,dimension(:,:), allocatable :: fun1,fun2
!     --------------------------------
      integer :: jj,kk
      integer :: iter,itermax
      integer :: nc1,nc2,nc3
      real*8  :: f1,f2,f3
      real*8  :: suma,sumaglob
!     --------------------------------
      real,dimension(2) :: tt
      real ::tcpu
      integer:: idmin,idhr,idsec
!     --------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         print*,' '
         write(*,'(t6,60a)'), '>>>>> Begin subroutine: TestSERIALnsmp3D'
         print*,' '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      IF (rang_topo.eq.0) THEN
         print*,' '
         write(*,'(t6,60a)'), '>>>>> Begin subroutine: TestMPInsmp3D'
         print*,' '
      ENDIF
#     endif
!     =============== END ================    
!     ====================================
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!      __________________________________
!     |                                  |
!     |             allocate             |
!     |__________________________________|

      allocate(fun1(N_CELL,NZ),fun2(N_CELL,NZ))

!      __________________________________
!     |                                  |
!     |        Initial function          |
!     |__________________________________|

      do k=1,NZ
         do i=1,N_CELL
            fun1(i,k) = k*dsin(xc(i)+yc(i))**2 
         enddo
      enddo

!*********************************************************************!
!                                                                     !
!                             Iterations                              !
!                                                                     !
!*********************************************************************!

      itermax = 3000
      iter = 0
 111  continue
      iter = iter + 1

!     ___________________________________
!     Global calculation (to compare)   

      suma = 0.0d0
      do k=1,NZ      
         do i=1,N_CELL0
            if (nbe(i).eq.0) then
	        nc1 = No_cp(i,1)
	        nc2 = No_cp(i,2)
	        nc3 = No_cp(i,3)
                f1  = fun1(nc1,k)
                f2  = fun1(nc2,k)
                f3  = fun1(nc3,k)
                fun2(i,k) = dsin(fun1(i,k) + f1 + f2 + f3)
            else
                fun2(i,k) = fun1(i,k)
            endif
            suma = suma + fun2(i,k)
         enddo
      enddo

!     ___________________________________
!     Communication   

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
          call SUM_parallel(suma,sumaglob)
          suma = sumaglob
          call communication3D(fun2)
#     endif
!     =============== END ================    
!     ====================================

!     ___________________________________
!     Update               
      do k=1,NZ      
         do i=1,N_CELL
            fun1(i,k) = fun2(i,k)  
         enddo
      enddo
!     ___________________________________
!     Criteria

      if (iter.lt.itermax) then
         goto 111
      else
         goto 112
      endif

112   continue
    
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      __________________________________
!     |                                  |
!     |         Display results          |
!     |__________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                        '
      print*,'     -------------------------------------------------- '
      print*,'                  Results Serial Version                '
      print*,'                   Iterations=',iterMAX
      print*,'     -------------------------------------------------- '
      print*,'                     '
      print*,'               sum = ',Suma
      print*,'                     '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                 Results in all procerssors             '
      print*,'                   Iterations=',iterMAX
      print*,'     ================================================== '
      print*,'                     '
      print*,'               sum = ',Suma
      print*,'                     '
      ENDIF
#     endif
!     =============== END ================    
!     ====================================

!      __________________________________
!     |                                  |
!     |    Displaying simulation time    |
!     |__________________________________|

917   format(t6,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
918   format(t6,'          user : ',f10.3,' sec')
919   format(t6,'        system : ',f10.3,' sec') 

!     -------------------------------------------------------
!     Calculating the simulation time     
      tcpu  = etime(tt)
      idhr  = tcpu/3600
      idmin = tcpu/60-idhr*60
      idsec = tcpu-(idhr*3600+idmin*60)

!     -------------------------------------------------------
!     Displaying 
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                        '
      print*,'     -------------------------------------------------- '
      print*,'                      Simulation Time                   '
      print*,'     -------------------------------------------------- '
      print*,'                                                        '
      write(*,917), tcpu,idhr,idmin,idsec
      print*,'                                                        '
      print*,'     -------------------------------------------------- '
      print*,'                                                        '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                      Simulation Time                   '
      print*,'     ================================================== '
      print*,'                                                        '
      ENDIF 
      write(*,917), tcpu,idhr,idmin,idsec
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                                                        '
      ENDIF
#     endif
!     =============== END ================    
!     ====================================
!      __________________________________
!     |                                  |
!     |            Deallocate            |
!     |__________________________________|

      deallocate(fun1,fun2)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(*,'(t6,60a)'), '<<<<< End   subroutine: TestSERIALnsmp3D'
         print*,' '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      IF (rang_topo.eq.0) THEN
         write(*,'(t6,60a)'), '<<<<< End   subroutine: TestMPInsmp3D'
         print*,' '
      ENDIF
#     endif
!     =============== END ================    
!     ====================================
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	         

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END of TestCPUnsmp3D                      !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
