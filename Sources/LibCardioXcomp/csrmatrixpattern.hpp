#ifndef CSRMATRIXPATTERN_HPP
#define CSRMATRIXPATTERN_HPP

#include <vector>
#include <iostream>
#include <petscsys.h>

namespace cardioxcomp
{
  //! \class CSRMatrixPattern
  //! \brief Sparse matrix pattern in Compressed Sparse Row (CSR) storage
  class CSRMatrixPattern {
    
  public:
    
    //! \brief Index of first nonzero in a row
    //! \param[in] rowIndex index of a row
    //! \return index of first nonzero in row of index rowIndex
    inline PetscInt & rowPointer(PetscInt rowIndex) {
#ifdef NDEBUG
      return m_rowptr[rowIndex];
#else
      return m_rowptr.at(rowIndex);
#endif
    }
    
    //! \brief Index of first nonzero in a row
    //! \param[in] rowIndex index of a row
    //! \return index of first nonzero in row of index rowIndex
    inline PetscInt const & rowPointer(PetscInt rowIndex) const {
#ifdef NDEBUG
      return m_rowptr[rowIndex];
#else
      return m_rowptr.at(rowIndex);
#endif
    }
    
    //! \brief Number of rows
    //! \return number of rows
    inline size_t numRows() const {
      return m_rowptr.size()-1;
    }
    
    //! \brief Number of nonzeros
    //! \return number of nonzeros
    inline size_t numNonzeros() const {
      return m_colind.size();
    }
    
    //! \brief Number of nonzeros in a row
    //! \param[in] rowIndex index of a row
    //! \return number of nonzeros in row of index rowIndex
    inline PetscInt numNonzerosInRow(size_t rowIndex) const {
#ifdef NDEBUG
      return m_rowptr[rowIndex+1] - m_rowptr[rowIndex];
#else
      return m_rowptr.at(rowIndex+1) - m_rowptr.at(rowIndex);
#endif
    }
    
    //! \brief Column index of a nonzero
    //! \param[in] nonzeroIndex index of a nonzero
    inline PetscInt const & columnIndex(size_t nonzeroIndex) const {
#ifdef NDEBUG
      return m_colind[nonzeroIndex];
#else
      return m_colind.at(nonzeroIndex);
#endif
    }
    
    //! \brief Column index of a nonzero
    //! \param[in] nonzeroIndex index of a nonzero
    inline PetscInt & columnIndex(size_t nonzeroIndex) {
#ifdef NDEBUG
      return m_colind[nonzeroIndex];
#else
      return m_colind.at(nonzeroIndex);
#endif
    }
    
    //! \brief Resizing
    //! \param[in] numRows number of rows
    //! \param[in] numNonzeros number of nonzeros
    inline void resize(size_t numRows,size_t numNonzeros) {
      m_rowptr.resize(numRows+1,0);
      m_colind.resize(numNonzeros,0);
    }
    
    //! \brief clear
    inline void clear() {
      m_rowptr.clear();
      m_colind.clear();
    }
    
    //! \param[in] rowPointers This vector contains row pointers.
    //!                        Its size equals numRows+1: the number of rows (numRows) plus one.
    //!                        The numRows first entries contain indices of first nonzero in each row.
    //!                        The last entry contains the number of nonzeros.
    //
    //! \param[in] columnIndices This vector contains column indices of nonzeros.
    //!                          Is size equals the number of nonzeros.
    //!                          Each entry contain the column index of the corresponding nonzero.
    //! \brief Modify pattern's content
    void set(std::vector<PetscInt> const & rowPointers,std::vector<PetscInt> const & columnIndices) {
      m_rowptr = rowPointers;
      m_colind = columnIndices;
    }
    
    //! \brief Accessor to column indices of nonzeros.
    //! \return Vector of column indices of nonzeros
    std::vector<PetscInt> const & columnIndices() const {
      return m_colind;
    }
    friend std::ostream& operator <<  (std::ostream& f, const CSRMatrixPattern & pattern );
    void print(int rankProc) const;
  private:
    
    std::vector< std::vector<PetscInt> > m_sparseRows;
    
    //! \brief sparse row pointers (index of first nonzero in each row)
    std::vector<PetscInt> m_rowptr;
    
    //! \brief column indices of nonzeros
    std::vector<PetscInt> m_colind;
    
  };
  
}


#endif // CSRMATRIXPATTERN_HPP
