!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               ---  TEST MPI COMMUNICATION 3D ---                    !
!                      Miguel Angel Uh Zapata                         !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

       SUBROUTINE TestMPInsmp3D(No_cp,nbe)

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

      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real,dimension(2) :: tt
      real ::tcpu
      integer:: idmin,idhr,idsec
!     --------------------------------
      integer :: jj,kk
      integer :: iter,itermax
      integer :: nc1,nc2,nc3
      real*8  :: f1,f2,f3
      real*8  :: suma,sumaglob,SumaNoMPI,SumaMPI
!     --------------------------------
      real*8 ,dimension(:,:), allocatable :: fun1_global,fun2_global
      real*8 ,dimension(:,:), allocatable :: fun1,fun2,fun2old
      real*8, dimension(:),   allocatable :: xc_global,yc_global 
      real*8, dimension(:),   allocatable :: xc_local,yc_local
      integer :: s,ielem,jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (rang_topo.eq.0) THEN
         write(*,'(t6,60a)'), '>>>>> Begin subroutine: TestMPI3D'
      ENDIF
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!     ====================================
!     =====  START PARALLEL OPTION =======
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'    ___________________________________________________ '
      print*,'   |                                                   |'
      print*,'   |       TEST EXAMPLE MPI 3D UNSTRUCTURED GRID       |'
      print*,'   |___________________________________________________|'
      print*,'                                                        '
      ENDIF
!     =============== END ================    
!     ====================================
!      __________________________________
!     |                                  |
!     |             allocate             |
!     |__________________________________|

      allocate(fun1_global(N_CELL0Global,NZ),fun2_global(N_CELL0Global,NZ))
      allocate(fun1(N_CELL,NZ),fun2(N_CELL,NZ),fun2old(N_CELL,NZ))
      allocate(xc_global(N_CELL0global))
      allocate(yc_global(N_CELL0global)) 
      allocate(xc_local(N_CELL))
      allocate(yc_local(N_CELL))
!      __________________________________
!     |                                  |
!     |            cell-centers          |
!     |__________________________________|

!     ___________________________________
!     global cell-center
      do i=1,N_CELL0global
	 jv1=No_vp_global(i,1)
	 jv2=No_vp_global(i,2)
	 jv3=No_vp_global(i,3)              
	 xc_global(i)=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
	 yc_global(i)=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
      enddo
!     ___________________________________
!     local cell-center
      do i=1,N_CELL
         ielem = Index_global(i)
         xc_local(i) = xc_global(ielem)
         yc_local(i) = yc_global(ielem)
      enddo   

!      __________________________________
!     |                                  |
!     |        Initial function          |
!     |__________________________________|

      do k=1,NZ
         do i=1,N_CELL0Global
            fun1_global(i,k) = k*dsin(xc_global(i)+yc_global(i))**2 
         enddo
      enddo

      do k=1,NZ
         do i=1,N_CELL
            fun1(i,k) = k*dsin(xc_local(i)+yc_local(i))**2  
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3D(fun1)
#     endif	
!     =============== END ================    
!     ====================================

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
         do i=1,N_CELL0global
            if (nbe_global(i).eq.0) then
	        nc1 = No_cp_global(i,1)
	        nc2 = No_cp_global(i,2)
	        nc3 = No_cp_global(i,3)
                f1  = fun1_global(nc1,k)
                f2  = fun1_global(nc2,k)
                f3  = fun1_global(nc3,k)
                fun2_global(i,k) = dsin(fun1_global(i,k) + f1 + f2 + f3)
            else
                fun2_global(i,k) = fun1_global(i,k)
            endif
            suma = suma + fun2_global(i,k)
         enddo
      enddo
      sumaNoMPI = suma

!     ___________________________________
!     New function calculation 

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
         enddo
      enddo 

!     ___________________________________
!     Communication          

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3D(fun2)
#     endif	
!     =============== END ================    
!     ====================================
  
!     ___________________________________
!     Sum of all elements       

      suma = 0.0d0
      do k=1,NZ      
         do i=1,N_CELL0
            suma = suma + fun2(i,k)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
          call SUM_parallel(suma,sumaglob)
          sumaMPI = sumaglob
#     else
          sumaMPI = suma 
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

      do k=1,NZ
         do i=1,N_CELL0global
            fun1_global(i,k) = fun2_global(i,k)
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
!     =====  START PARALLEL OPTION =======
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                 Results in all procerssors             '
      print*,'                   Iterations=',iterMAX
      print*,'     ================================================== '
      ENDIF
      print*,'     PROCESSOR:',rang_topo
      print*,'        sum        = ',suma
      print*,'        sum MPI    = ',sumaMPI
      print*,'        sum No MPI = ',SumaNoMPI
      print*,'                     '
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
!     =====  START PARALLEL OPTION =======
      IF (rang_topo.eq.0) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                      Simulation Time                   '
      print*,'     ================================================== '
      print*,'                                                        '
      print*,'      PROCESSOR:',rang_topo
      write(*,917), tcpu,idhr,idmin,idsec
      write(*,918), tt(1)
      write(*,919), tt(2)   
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                                                        '
      ENDIF
!     =============== END ================    
!     ====================================
!      __________________________________
!     |                                  |
!     |            Deallocate            |
!     |__________________________________|

      deallocate(fun1,fun2,fun2old)
      deallocate(fun1_global,fun2_global)
      deallocate(xc_global,yc_global)
      deallocate(xc_local,yc_local)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (rang_topo.eq.0) THEN
         write(*,'(t6,60a)'), '<<<<< End   subroutine: TestMPI3D'
         print*,' '
      ENDIF
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END of TestMPI                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
