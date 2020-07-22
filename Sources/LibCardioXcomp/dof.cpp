#include "stdafx.h"
#include <set>

#include "dof.hpp"

using namespace std;

namespace cardioxcomp {
  
  
#undef VERBOSE_DEBUG
  
  Dof::Dof(const Mesh& mesh,int numComp):
  m_mesh(mesh),m_numComp(numComp),m_numDof(mesh.nv*numComp)
  {};
  
  Dof::~Dof(){
    if(m_initAO) AODestroy(&m_ao);
    if(m_initMappingNodes) ISLocalToGlobalMappingDestroy(&m_mappingNodes);
  }
  
  void Dof::initializePattern(int numProc, int rankProc) {
    //
    //================================
    // Naive partition of the dof
    //================================
    //
    PetscInt numDofProc_tmp = m_numDof/numProc;
    PetscInt beginIdDofLocal = rankProc*numDofProc_tmp;
    PetscInt numDofProc;
    PetscInt endIdDofLocal;
    if(rankProc == numProc-1) {
      numDofProc = m_numDof - beginIdDofLocal;
      endIdDofLocal = m_numDof;
    } else {
      numDofProc = numDofProc_tmp;
      endIdDofLocal = beginIdDofLocal + numDofProc;
    }
#ifdef VERBOSE_DEBUG
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n========  Naive partition on Processor #%d  ========\n",rankProc);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"numDofProc = %d\n",numDofProc);
    int jverbose=0;
    for(int i=beginIdDofLocal;i<endIdDofLocal;i++){
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d -> %d\n",jverbose,i);
      jverbose++;
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"........ End of Naive partition on Processor #%d  ........\n",rankProc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
    //
    vector < set<int> > connected_nodes;
    connected_nodes.resize(numDofProc);
    int iglob,jglob;
    for (int it=0; it<m_mesh.nt; it++) {
      for(int i=0;i<3;i++){
        for(int icomp=0;icomp<m_numComp;icomp++){
          iglob = m_mesh(it,i) * m_numComp + icomp;
          if(iglob >= beginIdDofLocal && iglob < endIdDofLocal) {
            for(int j=0;j<3;j++){
              for(int jcomp=0;jcomp<m_numComp;jcomp++){
                jglob = m_mesh(it,j) * m_numComp + jcomp;
                if (iglob != jglob) connected_nodes[iglob-beginIdDofLocal].insert(jglob); // Parmetis does not want the diagonal
              }
            }
          }
        }
      }
    }
    //Storage in CSR format
    size_t numNonZero=0;
    for (size_t idof = 0; idof < numDofProc; idof++) {
      numNonZero += connected_nodes[idof].size();
    }
    m_pattern.resize(numDofProc,numNonZero);
    PetscInt pos = 0;
    for (PetscInt idof = 0; idof < numDofProc; idof++) {
      m_pattern.rowPointer(idof+1) = m_pattern.rowPointer(idof) + connected_nodes[idof].size();
      pos = 0;
      for (set<int>::const_iterator iCon = connected_nodes[idof].begin(); iCon != connected_nodes[idof].end(); iCon++) {
        m_pattern.columnIndex(m_pattern.rowPointer(idof) + pos) = *iCon;
        pos++;
      }
    }
    
  }
  
  
  void Dof::buildPattern(int rankProc, const PetscInt* dofRepartition){
    m_pattern.clear();
    
    vector < set<int> > connected_nodes;
    connected_nodes.resize(m_numDof);
    int iglob,jglob;
    for (int it=0; it<m_mesh.nt; it++) {
      for(int i=0;i<3;i++){
        for(int icomp=0;icomp<m_numComp;icomp++){
          iglob = m_mesh(it,i) * m_numComp + icomp;
          if (dofRepartition[iglob] == rankProc) {
            for(int j=0;j<3;j++){
              for(int jcomp=0;jcomp<m_numComp;jcomp++){
                jglob = m_mesh(it,j) * m_numComp + jcomp;
                connected_nodes[iglob].insert(jglob);
              }
            }
          }
        }
      }
    }
    
    //Storage in CSR format
    PetscInt numNonZero=0;
    for (int idof = 0; idof < m_numDof; idof++) {
      numNonZero += connected_nodes[idof].size();
    }
    m_pattern.resize(m_numDofLocal,numNonZero);
    PetscInt pos = 0;
    PetscInt cptDof = 0;
    for (int idof = 0; idof < m_numDof; idof++) {
      if (dofRepartition[idof] == rankProc) {
        m_pattern.rowPointer(cptDof+1) = m_pattern.rowPointer(cptDof) + connected_nodes[idof].size();
        pos = 0;
        for (set<int>::const_iterator iCon = connected_nodes[idof].begin(); iCon != connected_nodes[idof].end(); iCon++) {
          m_pattern.columnIndex(m_pattern.rowPointer(cptDof) + pos) = *iCon;
          pos++;
        }
        cptDof++;
      }
    }
  }
  
  
  /*!
   \brief use ParMetis to partition dof (line of the matrix).
   The function call parmetis to compute the dof partition which is stored in m_dofPart.
   Next, the application ordering (AO) used to fo from the global felisce ordering to the global petsc ordering
   is created.
   Finaly, the element partition is created and is store in m_eltPart.
   \param[in] numProc Number of process.
   \param[in] rankProc Rank of the current process.
   \param[in] comm MPI communicator.
   */
  void Dof::partitionDof(idx_t numProc,int rankProc,MPI_Comm comm) {
    /*!
     Partitioning of the dof.
     */
    
    /*
     Declaration of the parameters for ParMetis.
     idx_t define the type for ParMetis int (don't use int, keep idx_t !)
     */
    
    // Input parameters
    idx_t   num_flag = 0;                        // C-style array index.
    idx_t   ncon = 1;                            // number of weight per vertices of the graph.
    idx_t   wgtFlag = 0;                         // No weight in the graph.
    idx_t   option[3];                           // verbose level and random seed for parmetis.
    
    idx_t* vertexdist = new idx_t[numProc+1];

    real_t  ubvec = static_cast<real_t>(1.05);   // Imbalance tolerance
    real_t* tpwgts= new real_t[numProc*ncon];    // See parmetis manual for description
    
    // Output parameters
    idx_t   edgecut;                             // Number of edge in the graph that are cut by the partitioning.
    size_t numRows = m_pattern.numRows();
    std::vector<idx_t> dofPartitionning(numRows, 0);  // The partition of the dof.
    
    /*
     Initialisation of the parameters.
     */
    
    // vertexdist
    // vertexdist[j] = sum(size in process i), for i=0,..,j-1.
    // initialize the first value (always 0)
    vertexdist[0] = 0;
    
    // Gather all the size in vertexdist, starting at index 1 (because vertexdist[0] is 0)
    
    MPI_Allgather(&numRows, 1, MPI_INT, &vertexdist[1], 1, MPI_INT, comm);
    
    // Compute the partial sum
    for (int i = 1; i < numProc; i ++)
      vertexdist[i + 1] += vertexdist[i];
    
    // option
    option[0] = 1;     // 0 = default values, 1 = user defined, to print all that ugly information.
    option[1] = 0;   // verbose level.
    option[2] = 100;   // random seed.
    
    // tpwgts
    for ( PetscInt i= 0; i < numProc*ncon; i++)
      tpwgts[i] = static_cast<real_t>(1.)/numProc;
    
    /*
     Parmetis call
     */
    
    MPI_Comm comm_bis = PETSC_COMM_WORLD;

#ifdef PETSC_3_6
    idx_t mpatternRowPointer = m_pattern.rowPointer(0);
    idx_t mpatternColumnIndex = m_pattern.columnIndex(0);
    
    ParMETIS_V3_PartKway(vertexdist, &mpatternRowPointer, &mpatternColumnIndex, NULL, NULL, &wgtFlag, &num_flag, &ncon, &numProc, tpwgts, &ubvec, option, &edgecut, &dofPartitionning[0], &comm_bis);
#else
    ParMETIS_V3_PartKway(vertexdist, &m_pattern.rowPointer(0), &m_pattern.columnIndex(0), NULL, NULL, &wgtFlag, &num_flag, &ncon, &numProc, tpwgts, &ubvec, option, &edgecut, dofPartitionning.data(), &comm_bis);
#endif
    
    delete [] tpwgts;
    
    /*
     Gather the partitioning to all process
     */
    m_dofPart.resize(m_numDof, 0);
    
    int* recvcount = new int[numProc];//64bits
    for(int i=0; i<numProc; i++)
      recvcount[i] = vertexdist[i+1] - vertexdist[i];
    
#ifdef PETSC_3_6
    MPI_Allgatherv(&dofPartitionning[0], dofPartitionning.size(), MPI_INT, &m_dofPart[0], recvcount, vertexdist, MPI_INT, comm);
#else
    MPI_Allgatherv(dofPartitionning.data(), dofPartitionning.size(), MPI_INT, m_dofPart.data(), recvcount, vertexdist, MPI_INT, comm);    
#endif

    delete [] recvcount;
    
    // Count the number of local dofs
    m_numDofLocal = 0;
    for (PetscInt i = 0; i < m_numDof; i++) {
      if (m_dofPart[i] == rankProc) m_numDofLocal++;
    }
    
#ifdef VERBOSE_DEBUG
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n======== Print Parmetis partition on Processor #%d  ========\n",rankProc);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"numDofLocal = %d\n",m_numDofLocal);
    int jtmp=0;
    for (PetscInt i = 0; i < m_numDof; i++){
      if (m_dofPart[i] == rankProc){
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d -> %d\n",jtmp,i);
        jtmp++;
      }
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"......... End of Parmetis partition on Processor #%d ..........\n",rankProc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
    /*
     Build the matrix pattern
     */

#ifdef PETSC_3_6
    vector<platformInt> mDofPart;
    for(unsigned int il=0 ; il<m_dofPart.size() ; il++){
      mDofPart.push_back(m_dofPart[il]);
    }
    buildPattern(rankProc, mDofPart.data());
#else
    buildPattern(rankProc, m_dofPart.data());
#endif
    
    /*!
     Global ordering: AO is an application between global mesh ordering and global ordering of Petsc
     AO: globalordering to globalordering
     the action is:
     - Take the global index associated to the i-th proc
     - Create an application to a new global index in order to have contiguous global indexes on the processor
     
     both global ordering are on input (loc2glob and pordering)
     AOCreateBasic(comm, m_numDofLocal, loc2Glob, pordering, &m_ao);
     
     Example:
     Proc 0:
     glob2loc = [0 1 2 3 8]
     pordering =[0 1 2 3 4]
     
     Proc 1:
     glob2loc = [4 5 6 7]
     pordering =[5 6 7 8]
     */
    
    platformInt rstart;
    PetscInt cptDof = 0;
    
    PetscInt* loc2Glob = new PetscInt[m_numDofLocal];
    for (PetscInt i = 0; i < m_numDofLocal; i++) loc2Glob[i] = 0;
    
    for (PetscInt i = 0; i < m_numDof; i++) {
      if (m_dofPart[i] == rankProc) {
        // cptDof (local) -> i (global)
        loc2Glob[cptDof] = i;
        cptDof++;
      }
    }
    
    // inclusive sum (rank_0 = numdoflocal_0, rank_i = sum_j=0^rank_i numdoflocal_i)
    MPI_Scan(&m_numDofLocal, &rstart, 1, MPIU_INT, MPI_SUM, comm);
    rstart -= m_numDofLocal;

#ifdef PETSC_3_6
    int value;
    MPI_Scan(&m_numDofLocal, &value, 1, MPIU_INT, MPI_SUM, comm);
    value -= m_numDofLocal;
#endif
    
    PetscInt *pordering = new PetscInt[m_numDofLocal];
    for (PetscInt i= 0; i < m_numDofLocal; i++) {
      // pordering: so global indexes are contiguous on a processor; pordering: global->global
      pordering[i] = rstart + i;
    }
    
    AOCreateBasic(comm, m_numDofLocal, loc2Glob, pordering, &m_ao);
#ifdef VERBOSE_DEBUG
    AOView(m_ao,PETSC_VIEWER_STDOUT_WORLD);
#endif
    m_initAO = true;
    
    delete [] loc2Glob;
    delete [] pordering;
    
    /*!
     Compute the element repartition.
     An element is owned by a process k if most of its dof are on process k.
     In case of a draw, the element belongs to the process with the smallest id.
     */
    
    PetscInt max = 0;                                  // used to compute indice_max.
    PetscInt indice_max = 0;                           // id of the process that owns the element.
    PetscInt idDof;                                    // Local number of a dof.
    PetscInt *numDofPerProcs = new PetscInt[numProc];    // Number of dofs per process. //per element to be computed for each element and then cleared
    
    for ( int i = 0; i < numProc; i++) numDofPerProcs[i] = 0;
    m_eltPart.resize(m_mesh.nt, 0);
    for(int it=0;it<m_mesh.nt;it++){
      for(int i=0;i<3;i++){
        for(int icomp=0;icomp<m_numComp;icomp++){
          idDof = m_mesh(it,i) * m_numComp + icomp;
          numDofPerProcs[m_dofPart[idDof]]++;
        }
      }
      // Find the rank of the process that owns the most dof in the current element.
      for (int j = 0; j < numProc; j++) {
        if (max < numDofPerProcs[j]) {
          max = numDofPerProcs[j];
          indice_max = j;
        }
      }
      m_eltPart[it] = indice_max;
      // Clear local variables
      max = 0;
      indice_max = 0;
      for ( int j = 0; j < numProc; j++)
        numDofPerProcs[j] = 0;
#ifdef VERBOSE_DEBUG
      PetscPrintf(PETSC_COMM_WORLD,"Element %d belongs to proc #%d\n", it, m_eltPart[it]);
#endif
    }
  }
  
  void Dof::allocateMatrixVec(int numProc, int rankProc,Mat* matrix,int numMat,Vec* vec,int numVec){
    std::cout << "Allocate matrices and vectors on proc #" << rankProc << std::endl;
    // Create the local to global mapping for support dofs.
    PetscInt *loc2Glob = new PetscInt[m_numDofLocal];
    PetscInt cptDof= 0;
    for ( PetscInt iLDof = 0; iLDof < m_numDofLocal; iLDof++)
      loc2Glob[iLDof] = 0;
    
    for ( PetscInt iDof = 0; iDof < m_numDof; iDof++) {
      if ( m_dofPart[iDof] == rankProc ) {
        loc2Glob[cptDof] = iDof;
        cptDof++;
      }
    }
    
    // global to local mapping on nodes
#ifdef PETSC_3_6
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, m_numDofLocal,1,loc2Glob, PETSC_COPY_VALUES, &m_mappingNodes);
#else
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, m_numDofLocal, loc2Glob, PETSC_COPY_VALUES, &m_mappingNodes);
#endif
    // ISLocalToGlobalMappingView( m_mappingNodes, PETSC_VIEWER_STDOUT_WORLD);
    delete [] loc2Glob;
    m_initMappingNodes = true;
    if ( numProc > 1 ) {
      //
      // Multi processors
      //
      PetscInt* iCSR = new PetscInt[m_pattern.numRows()+1];
      PetscInt* jCSR = new PetscInt[m_pattern.numNonzeros()];
      iCSR[0] = 0;
      set<PetscInt> jCSRSorted;
      PetscInt connect;
      PetscInt globValue;
      for ( PetscInt ii = 0; ii < m_numDofLocal; ii++) {
        connect = m_pattern.numNonzerosInRow(ii);
        iCSR[ii+1] = iCSR[ii] + connect;
        for ( PetscInt jj = 0; jj < connect; jj++) {
          globValue = m_pattern.columnIndex( m_pattern.rowPointer(ii) + jj );
          AOApplicationToPetsc(m_ao, 1, &globValue);
          jCSRSorted.insert( globValue );
        }
        PetscInt kk = 0;
        for ( set<PetscInt>::const_iterator  it_jCSR = jCSRSorted.begin(); it_jCSR != jCSRSorted.end(); it_jCSR++) {
          jCSR[ iCSR[ii] + kk] = *it_jCSR;
          kk++;
        }
        jCSRSorted.clear();
      }
      for(int i=0;i<numMat;i++){
        MatCreateAIJ(PETSC_COMM_WORLD, m_numDofLocal, m_numDofLocal, m_numDof, m_numDof,0, PETSC_NULL, 0, PETSC_NULL, &matrix[i]);
        MatSetFromOptions(matrix[i]);
        MatMPIAIJSetPreallocationCSR(matrix[i], iCSR, jCSR, NULL);
        MatSetOption(matrix[i], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      }
      
      for(int i=0;i<numVec;i++){
        VecCreateMPI(PETSC_COMM_WORLD, m_numDofLocal, m_numDof, &vec[i]);
      }
      delete [] iCSR;
      delete [] jCSR;
    } else {
      //
      // One processor
      //
      size_t numRows = m_pattern.numRows();
      PetscInt* nnz = new PetscInt[ numRows ];
      for ( PetscInt idof = 0; idof < static_cast<PetscInt>(numRows); idof++) {
        nnz[idof] = m_pattern.numNonzerosInRow(idof);
      }
      
      for(int i=0;i<numMat;i++){
        MatCreateSeqAIJ(PETSC_COMM_WORLD, numRows, numRows, 0, nnz, &matrix[i]);
        MatSetFromOptions(matrix[i]);
      }
      for(int i=0;i<numMat;i++){
        MatSetOption(matrix[i], MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      }
      
      for(int i=0;i<numVec;i++){
        VecCreate(PETSC_COMM_WORLD, &vec[i]);
        VecSetSizes(vec[i], numRows, numRows );
        VecSetFromOptions(vec[i]);
      }      
      delete [] nnz;
    }
  std::cout << "End allocate matrices and vectors on proc #" << rankProc << std::endl;
  }
}
