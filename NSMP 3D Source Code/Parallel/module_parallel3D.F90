!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             MODULE PARALLEL 3D PROBLEM (UNSTRUCTURED MESH)          !
!                      Miguel Angel Uh Zapata                         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:     1.1)  initialisation_mpi                      !
!                       1.2)  parallel_neighbors                      !
!                       1.3)  parallel_topology                       !
!                       1.4)  parallel_index                          !
!                       1.5)  parallel_shareindex                     !
!                       1.6)  finalisation_mpi                        !
!                                                                     !
!---------------------------------------------------------------------!

      MODULE parallel

      implicit none
#     include "mpif.h"
#     include "common.mpf"
!     ------------------------------------
      integer :: code,rang_topo
      integer :: comm3D,comm2
      integer :: Nprocs,Nindex,Nedges
      integer :: TotalSubdomains
      integer :: NeighborNumber
      integer :: NextraMax,NshareMax
!     ------------------------------------
      integer,dimension(:), allocatable :: edges
      integer,dimension(:), allocatable :: indexval
      integer,dimension(:), allocatable :: Neighbors
!     ------------------------------
      integer,dimension(:), allocatable :: type_send
      integer,dimension(:), allocatable :: type_recv
      integer,dimension(:), allocatable :: type_sendVEC
      integer,dimension(:), allocatable :: type_recvVEC
      integer :: type_bloc1
!     ------------------------------------
      integer,dimension(:,:), allocatable :: CCdom
      integer,dimension(:,:), allocatable :: DimCCdom
      integer,dimension(:,:), allocatable :: ExtraCCdom
      integer,dimension(:,:), allocatable :: DimExtraCCdom
      integer,dimension(:,:), allocatable :: VVdom
      integer,dimension(:,:), allocatable :: DimVVdom
!     ------------------------------
      integer,dimension(:,:), allocatable :: IniExtraFrom
      integer,dimension(:,:), allocatable :: FinExtraFrom
      integer,dimension(:,:), allocatable :: DimExtraFrom
!     ------------------------------
      integer,dimension(:),   allocatable :: Index_global
      integer,dimension(:),   allocatable :: Index_globalv
      integer,dimension(:),   allocatable :: IniExtraIndex_local
      integer,dimension(:,:), allocatable :: ShareIndex_local
      integer,dimension(:),   allocatable :: ShareIndex_dim
      real*8, dimension(:),   allocatable :: SharePhi
      real*8, dimension(:,:), allocatable :: SharePhiNew
!     ------------------------------
      integer,dimension(:),   allocatable :: wb_aux,qb_aux,hb_aux,sp_aux
!     ------------------------------
      integer,dimension(:),   allocatable :: TagParallel
!     ------------------------------
      integer,dimension(:,:), allocatable :: No_vp_global
      integer,dimension(:,:), allocatable :: No_cp_global
      integer,dimension(:),   allocatable :: nbe_global
      integer,dimension(:),   allocatable :: No_wb_global
      integer,dimension(:),   allocatable :: No_qb_global
      integer,dimension(:),   allocatable :: No_hb_global
      integer,dimension(:),   allocatable :: No_sp_global
      real*8, dimension(:),   allocatable :: xv_global
      real*8, dimension(:),   allocatable :: yv_global
      real*8, dimension(:),   allocatable :: zbv_global

      CONTAINS
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.1  |                 Initialization MPI                     |
!     |______|________________________________________________________|

      SUBROUTINE initialisation_mpi

      !--------------------------------------------------------------!   
      !                                                              !
      !   ---> code      : name of the MPI program                   !
      !   <--- Nprocs    : total number of processors                !
      !   <--- rang_topo : rank (name) of each processor             !
      !                                                              !
      !--------------------------------------------------------------!
 
!     __________________________
!     Initialisation de MPI
      call MPI_INIT(code)
!     __________________________
!     Number of processors 
      call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,code) 
!     __________________________
!     Rang of each processor
      call MPI_COMM_RANK(MPI_COMM_WORLD,rang_topo,code)

!     __________________________
!     Display global parameters
      IF (rang_topo.eq.0) THEN   
      print*,'                                                         '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'=========================================================' 
      print*,'           __________________________________            '
      print*,'          |                                  |           '
      print*,'          |  ***  MPI PROGRAM: NSMP3D   ***  |           '
      print*,'          |__________________________________|           '
      print*,'                                                         '
      print*,'========================================================='
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                         '
      print*,'   ====================================================  '
      print*,'                   MPI: INITIALIZATION                   '
      print*,'   ====================================================  '
      print*,'                                                         '
      print*,'       NUMBER OF MPI PROCESSORS: ', Nprocs
      print*,'                                                         '
      print*,'   ____________________________________________________  '
      print*,'                                                         '
      print*,'       GLOBAL DOMAIN:                                    ' 
      print*,'          Global number of elements  =',N_CELL0Global 
      print*,'          Global number of vertices  =',N_VERTGlobal
      print*,'          Global number of NZ points =',NZglobal-1
      print*,'                                                         '
      ENDIF

      END SUBROUTINE initialisation_mpi

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.2  |               Processors topology                      |
!     |______|________________________________________________________|

      SUBROUTINE parallel_topology

      !---------------------------------------------------------------!   
      !                                                               !
      !  Division of the domain in subdomains by a graph topology.    !
      !                                                               !
      ! ---> MPI_COMM_WORLD : the communicator group we are using.    !
      ! ---> nbprocs: the number of processors.                       !
      ! ---> index  : array of integers describing node degrees       !
      ! ---> edges  : array of integers describing graph edges        !
      ! --->   0    : if we don't want to order processes in the group!
      ! <--- comm3D : the communicator which represents the graph     !
      !                                                               !
      !---------------------------------------------------------------!

      call MPI_GRAPH_CREATE(MPI_COMM_WORLD,Nprocs,indexval,edges,0,comm3D,code)

!     __________________________
!     Display topology
      IF (rang_topo.eq.1) THEN  
      print*,'   ____________________________________________________  '
      print*,'                                                         '
      print*,'       MPI GRAPH TOPOLOGY:                               '
      print*,'          Number of indexes =',Nindex
      print*,'          Number of edges   =',Nedges
     !print*,'          Index values      =',indexval(0:Nindex-1)
     !print*,'          Edge values       =',edges(0:Nedges-1)
      print*,'                                                        '
      ENDIF

      END SUBROUTINE parallel_topology

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.3  |                    Neighbors                           |
!     |______|________________________________________________________|

       SUBROUTINE parallel_neighbors 

      !---------------------------------------------------------------!   
      !                                                               !
      !  Calculation of the neighbors:                                !
      !  <--- neighbourNumber : is number of neighbors of the "s" proc!
      !  <--- neighbour       : is the array neighbors of "s" proc.   !
      !                                                               !
      !---------------------------------------------------------------!   

!     _________________________________________________________________
!     Neighbors calculation
!     --------------------------------
!     Number of neighborns 
      call MPI_Graph_neighbors_count(comm3D,rang_topo,NeighborNumber,code)
!     --------------------------------
!     Array of neighborns 
      allocate(Neighbors(NeighborNumber))
      call MPI_Graph_neighbors(comm3D,rang_topo,NeighborNumber,Neighbors)

!     _________________________________________________________________
!     Display

      IF (rang_topo.eq.1) THEN  
      print*,'   ____________________________________________________  '
      print*,'                                                         '
      print*,'       NEIGHBORS:                                        '
      ENDIF
      print*,'          Processor:',rang_topo
      print*,'          Number of neighbors =',NeighborNumber
      print*,'          Processor neighbors :',Neighbors(1:NeighborNumber)
      print*,'          ----------------------------------'
      print*,'                                                          '

      END SUBROUTINE parallel_neighbors


!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.4  |           Assignation of global index                  |
!     |______|________________________________________________________|

      SUBROUTINE parallel_index

!---------------------------------------------------------------------!
!                                                                     !
!     Global index values corresponding to local index values for     !
!     each subdomain.                                                 !
!                                                                     !
!     IMPORTANT: It is really important to remark the order of the    !
!                overlaping cell-centers are located. I will include  !
!                the ghost points and later the extra points.         !
!                                                                     !
!---------------------------------------------------------------------!

      integer :: Iini,Ifin,i0,s
      integer :: k,ii,nv

      s = rang_topo + 1

      allocate(Index_globalv(N_VERT))
      allocate(Index_global(N_CELL))
!     ________________________________________________________
!     Global index of vertex points in each subdomain     

      nv = 0
      Iini = DimVVdom(s,2)
      Ifin = DimVVdom(s,3)
      do i=Iini,Ifin
         nv = nv + 1
         Index_globalv(nv) = VVdom(i,1)
      enddo
!      ________________________________________________________
!     Global index of cell-center points in each subdomain 

!     -----------------
!     Interior 
      i0 = 0
      Iini = DimCCdom(s,2)
      Ifin = DimCCdom(s,3)
      do i=Iini,Ifin
         i0 = i0 + 1
         Index_global(i0) = CCdom(i,1)
      enddo
!     -----------------
!     Ghost (No global index needed, temporal = -10) 
      do i=1,N_CELLghost
         i0 = i0 + 1
         Index_global(i0) = -10
      enddo      
!     -----------------
!     extra overlaping  
      Iini = DimExtraCCdom(s,2)
      Ifin = DimExtraCCdom(s,3)
      do i=Iini,Ifin
         i0 = i0 + 1
         Index_global(i0) = ExtraCCdom(i,3)
      enddo

!     ________________________________________________________
!     Free memory CCdom global

      deallocate(CCdom)
      deallocate(DimCCdom)

      END SUBROUTINE parallel_index

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.5  |              Index to communication                    |
!     |______|________________________________________________________|

      SUBROUTINE parallel_shareindex

      !---------------------------------------------------------------!   
      !                                                               !
      !       Assign the vertex and cell-centers to each domain       !
      !                                                               !
      !   <--- IniExtraIndex_local: Initial local index to received   !
      !   <--- ShareIndex_local   : Index send to neighbors           !
      !   <--- ShareIndex_dim     : Number of index send              !  
      !   <--- SharePhi           : Just allocate this vector         !
      !                                                               !
      !---------------------------------------------------------------!  

      integer :: Iini,Ifin,i0,ielem,jelem,ivert,s,sm,suma
      integer :: k,ii
      integer :: nv1,nv2,nv3

      s = rang_topo + 1

!     ________________________________________________________
!     Initial index of extra cells by neighborn 

      allocate(IniExtraIndex_local(NeighborNumber))

      s  = rang_topo + 1
      suma = N_CELL0 + N_CELLghost + 1
      do i=1,NeighborNumber
          sm = Neighbors(i) + 1
          IniExtraIndex_local(i) = suma
          suma = suma + DimExtraFrom(s,sm)
      enddo
!     ________________________________________________________
!     Share index to send

      NextraMax = 0
      do s=1,TotalSubdomains
         NextraMax = max(NextraMax,DimExtraCCdom(s,1))
      enddo

      allocate(ShareIndex_local(NextraMax,NeighborNumber))
      allocate(ShareIndex_dim(NeighborNumber))

      s = rang_topo + 1
      do k=1,NeighborNumber
         sm = Neighbors(k) + 1
         Iini = IniExtraFrom(sm,s)
         Ifin = FinExtraFrom(sm,s)
         i0 = 0
         do j=Iini,Ifin
            i0 = i0 + 1
            jelem = ExtraCCdom(j,3)
            do ii=1,N_CELL0
               ielem = Index_global(ii)
               if (ielem.eq.jelem) then
                   ShareIndex_local(i0,k) = ii
               endif
            enddo
         enddo
         ShareIndex_dim(k) = i0
      enddo

!     ________________________________________________________
!     Free memory ExtraCCdom global

      deallocate(ExtraCCdom)
      deallocate(DimExtraCCdom)

      END SUBROUTINE parallel_shareindex

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.6  |                 Finalization MPI                       |
!     |______|________________________________________________________|

      SUBROUTINE finalisation_mpi

      !---------------------------------------------------------------!   
      !                                                               !
      !   Finalization of definitions and the MPI version             !
      !                                                               !
      !---------------------------------------------------------------! 

!     ________________________________________________________
!     Free memory
!     __________________
!     Index and edges 
      deallocate(indexval,edges)
      deallocate(Neighbors)
!     __________________
!     Global variables
!     ---------------------
      deallocate(VVdom) 
      deallocate(DimVVdom)
!     ---------------------
      deallocate(IniExtraFrom)
      deallocate(FinExtraFrom)
      deallocate(DimExtraFrom)
!     ---------------------
      deallocate(No_vp_global,  &
                 No_cp_global,  & 
                 nbe_global,    &      
                 No_wb_global,  & 
                 No_qb_global,  &
                 No_hb_global,  &
                 No_sp_global,  &              
                 xv_global,     &
                 yv_global,     &
                 zbv_global)
!     __________________
!     Local variables
      deallocate(Index_global)
      deallocate(Index_globalv)
!     ---------------------
      deallocate(IniExtraIndex_local)
      deallocate(ShareIndex_local)
      deallocate(ShareIndex_dim)
      deallocate(SharePhi)
      deallocate(SharePhiNew)
!     __________________
!     Communication variables
      deallocate(type_recv)
      deallocate(type_send)
      deallocate(type_recvVEC)
      deallocate(type_sendVEC)
!     ________________________________________________________	
!     Free block vertex type 
      call MPI_TYPE_FREE(type_bloc1,code)

!     ________________________________________________________	
!     Free code
      call MPI_FINALIZE(code)


      END SUBROUTINE finalisation_mpi

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   MODULE PARALLEL (COMMUNICATION)                   !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:    3.1)  parallel_type                            !
!                      3.1)  communication2D                          !
!                      3.2)  communication3D                          !
!                      3.4)  communication3Dnew                       !
!                                                                     !
!---------------------------------------------------------------------!

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.1  |    Definition of the type of communication vectors     |
!     |______|________________________________________________________|

      SUBROUTINE parallel_type

      !---------------------------------------------------------------!   
      !                                                               !
      ! Specifying the vector type:                                   !
      !                                                               !
      !   <--- type_send(i)   : vector of extra elements to send to   !
      !                         each neighborn i (continuous type).   !
      !   <--- type_recv(i)   : vector of extra elements to receive to!
      !                         each neighborn i (continuous type).   !
      !   <--- type_sendVEC(i): vector of extra elements to send to   !
      !                         each neighborn i (vector type).       !
      !   <--- type_recvVEC(i): vector of extra elements to receive to!
      !                         each neighborn i (vector type).       !
      !                                                               !
      !---------------------------------------------------------------! 

      integer :: number_send,number_recv
      integer :: i,s,sm

      allocate(type_send(NeighborNumber))
      allocate(type_recv(NeighborNumber))
      allocate(type_sendVEC(NeighborNumber))
      allocate(type_recvVEC(NeighborNumber))

!     ________________________________________
!     Continuous type: size of sent array
      do i=1,NeighborNumber     
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_send = DimExtraFrom(sm,s)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_send,&
              MPI_DOUBLE_PRECISION,type_send(i),code)
         call MPI_TYPE_COMMIT(type_send(i),code)
!        --------------------
!        Vector type
         call MPI_TYPE_VECTOR(NZ,number_send,N_CELL,&
              MPI_DOUBLE_PRECISION,type_sendVEC(i),code)
         call MPI_TYPE_COMMIT(type_sendVEC(i),code)
      enddo

!     ________________________________________
!     Continuous type: size of received array
      do i=1,NeighborNumber 
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_recv = DimExtraFrom(s,sm)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_recv,&
              MPI_DOUBLE_PRECISION,type_recv(i),code)
         call MPI_TYPE_COMMIT(type_recv(i),code)
!        --------------------
!        Vector type
         call MPI_TYPE_VECTOR(NZ,number_recv,N_CELL,&
              MPI_DOUBLE_PRECISION,type_recvVEC(i),code)
         call MPI_TYPE_COMMIT(type_recvVEC(i),code)
      enddo

!     ________________________________________
!     Vertex global block type      
      call MPI_TYPE_VECTOR(NZglobal-1,N_VERTglobal,N_VERTglobal,&
           MPI_DOUBLE_PRECISION,type_bloc1,code)
      call MPI_TYPE_COMMIT(type_bloc1,code)

!     ________________________________________
!     Allocate auxiliar variable to communicate 

      NshareMax = 0
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,ShareIndex_dim(i))
      enddo

      allocate(SharePhi(NshareMax))
      allocate(SharePhiNew(N_CELL,NZ))


      END SUBROUTINE parallel_type

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.2  |        Communication continuos type 2D                 |
!     |______|________________________________________________________|

       SUBROUTINE communication2D(phi)

      !---------------------------------------------------------------!   
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable ar with one entry phi(N_CELL).                    !
      !                                                               !
      !---------------------------------------------------------------! 

      real*8,dimension(:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,elem,Iini 

!     ______________________________ 
!     Send information
      do i=1,NeighborNumber
!        -----------
         do j=1,ShareIndex_dim(i)
            elem = ShareIndex_local(j,i)
            SharePhi(j) = phi(elem)
         enddo
!        -----------
         call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                       etiquette,comm3D,code)
      enddo

      call MPI_Barrier(comm3D,code)
!     ______________________________ 
!     Receive information
      do i=1,NeighborNumber
         Iini = IniExtraIndex_local(i)
         call MPI_RECV(phi(Iini),1,type_recv(i),Neighbors(i),&
                       etiquette,comm3D,statut,code)
      enddo

      END SUBROUTINE communication2D

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.3  |         Communication 3D: vector type                  |
!     |______|________________________________________________________|

       SUBROUTINE communication3Dtype1(phi)

      !---------------------------------------------------------------!   
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !                                                               !
      !---------------------------------------------------------------! 

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini 

!     ______________________________ 
!     Send information
      do i=1,NeighborNumber
!        -----------
         do j=1,ShareIndex_dim(i)
            elem = ShareIndex_local(j,i)
            do k=1,NZ
               SharePhiNew(j,k) = phi(elem,k)
            enddo
         enddo
!        -----------
         call MPI_SEND(SharePhiNew(1,1),1,type_sendVEC(i),Neighbors(i),&
                       etiquette,comm3D,code)
      enddo
      call MPI_Barrier(comm3D,code)
!     ______________________________ 
!     Receive information
      do i=1,NeighborNumber
         Iini = IniExtraIndex_local(i)
         call MPI_RECV(phi(Iini,1),1,type_recvVEC(i),Neighbors(i),&
                       etiquette,comm3D,statut,code)
      enddo


      END SUBROUTINE communication3Dtype1

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.4  |        Communication 3D: NZ continuos (type 2)         |
!     |______|________________________________________________________|

       SUBROUTINE communication3D(phi)

      !---------------------------------------------------------------!   
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !                                                               !
      !---------------------------------------------------------------! 

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini 

!     ______________________________ 
!     Send information
      DO k=1,NZ
         do i=1,NeighborNumber
!           -----------
            do j=1,ShareIndex_dim(i)
               elem = ShareIndex_local(j,i)
               SharePhi(j) = phi(elem,k)
            enddo
!           -----------
            call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                          etiquette,comm3D,code)
         enddo
      ENDDO

      call MPI_Barrier(comm3D,code)
!     ______________________________ 
!     Receive information
      DO k=1,NZ
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,k),1,type_recv(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
      ENDDO

      END SUBROUTINE communication3D


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           MODULE PARALLEL (ESPECIAL FUNCTIONS: min,max,sum)         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:     3.1)  MIN_parallel                            !
!                       3.2)  MAX_parallel                            !
!                       3.3)  SUM_parallel                            !
!                                                                     !
!---------------------------------------------------------------------!

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.1  |              MIN value all processors                  |
!     |______|________________________________________________________|

       SUBROUTINE MIN_parallel(value,valuemin)

      !---------------------------------------------------------------!   
      !                                                               !
      !    Return the minimum value among all subdomains:             !
      !    <--- valuemin = min{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------! 

       real*8 :: value,valuemin

       call MPI_ALLREDUCE(value,valuemin,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
                          comm3D,code)

       END SUBROUTINE MIN_parallel

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.2  |               MAX value all processors                 |
!     |______|________________________________________________________|

       SUBROUTINE MAX_parallel(value,valuemax)

      !---------------------------------------------------------------!   
      !                                                               !
      !    Return the maximum value among all subdomains:             !
      !    <--- valuemax = max{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------! 

       real*8 :: value,valuemax

       call MPI_ALLREDUCE(value,valuemax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                          comm3D,code)

       END SUBROUTINE MAX_parallel

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 3.3  |                 SUM value all processors               |
!     |______|________________________________________________________|

       SUBROUTINE SUM_parallel(value,valuesum)

      !---------------------------------------------------------------!   
      !                                                               !
      !    Return the addiction of a value among all subdomains:      !
      !    <--- valuesum = sum{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------! 

       real*8 :: value,valuesum

       call MPI_ALLREDUCE(value,valuesum,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                          comm3D,code)

       END SUBROUTINE SUM_parallel

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   MODULE PARALLEL (GLOBAL MATRIX)                   !
!                              Jul 2014                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:    4.1)  matgloV                                  !
!                      4.2)  matgloC                                  !
!                                                                     !
!---------------------------------------------------------------------!

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.1  |           Building the global matrix (vertices)        |
!     |______|________________________________________________________|

      SUBROUTINE matgloV(phiv,phiv_global)

      !---------------------------------------------------------------!   
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !                                                               !
      !---------------------------------------------------------------! 

      real*8, dimension(:,:) :: phiv,phiv_global
      integer,parameter :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,vert,ii,Iini,Ifin
      real*8,dimension(:,:),allocatable :: arglo
      

      allocate(arglo(N_VERTglobal,NZglobal-1))
!      ___________________________________________
!     |                                           |
!     |                ONE PROCESSOR              |
!     |___________________________________________|

      IF (Nprocs.eq.1) THEN
         do k=1,NZ-1
            do nv=1,N_VERT
               vert = Index_globalv(nv)
               phiv_global(vert,k) = phiv(nv,k)
            enddo
         enddo
!      ___________________________________________
!     |                                           |
!     |          MORE THAN ONE PROCESSOR          |
!     |___________________________________________|

      ELSE
!       ___________________________________	
!       SENDIND INFORMATION 
        if (rang_topo.eq.0) then
           do k=1,NZ-1
              do nv=1,N_VERT
                 vert = Index_globalv(nv)
                 phiv_global(vert,k) = phiv(nv,k)
              enddo
           enddo
        else
           do k=1,NZ-1
              do nv=1,N_VERT
                 arglo(nv,k) = phiv(nv,k)
              enddo
           enddo
	   call MPI_SEND(arglo(1,1),1,type_bloc1,0,&
                         etiquette,comm3D,code)
        endif
!       ___________________________________	
!       RECEIVING INFORMATION 
        if (rang_topo.eq.0) then
	    do s=1,Nprocs-1
	       call MPI_RECV(arglo(1,1),1,type_bloc1,s,&
     		             etiquette,comm3D,statut,code)
               Iini = DimVVdom(s+1,2)
               Ifin = DimVVdom(s+1,3)
               nv = 0
               do i=Iini,Ifin
                  nv = nv + 1
                  vert = VVdom(i,1)
                  do k=1,NZglobal-1
                     phiv_global(vert,k) = arglo(nv,k)
                  enddo
               enddo
	    enddo
        endif
      ENDIF

      deallocate(arglo)

      END SUBROUTINE matgloV

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 4.2  |   Building the global matrix (cell-centers) NO YET !!! |
!     |______|________________________________________________________|

      SUBROUTINE matgloC(phi,phi_global)

      !---------------------------------------------------------------!   
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phi_global(N_CELL0global)                   !
      !    using all the subdomain matrices phi(N_CELL0).             !
      !                                                               !
      !---------------------------------------------------------------! 

      real*8, dimension(:,:) :: phi,phi_global
      integer:: i,k,elem


      END SUBROUTINE matgloC

      END MODULE parallel

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END MODULE PARALLEL                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
