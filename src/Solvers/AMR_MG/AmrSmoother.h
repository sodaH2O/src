#ifndef AMR_SMOOTHER_H
#define AMR_SMOOTHER_H

#include <string>

#include "AmrMultiGridDefs.h"

#include "Ifpack2_Factory.hpp"

/// Interface to Ifpack2 smoothers of the Trilinos package
class AmrSmoother {
    
public:
    typedef amr::global_ordinal_t   go_t;
    typedef amr::local_ordinal_t    lo_t;
    typedef amr::scalar_t           scalar_t;
    typedef amr::node_t             node_t;
    typedef amr::matrix_t           matrix_t;
    typedef amr::vector_t           vector_t;
    typedef Ifpack2::Preconditioner<scalar_t,
                                    lo_t,
                                    go_t,
                                    node_t
            > preconditioner_t;
    
    /// All supported Ifpack2 smoothers
    enum Smoother {
        GAUSS_SEIDEL = 0,
        SGS,    // symmetric Gauss-Seidel
        JACOBI //,
//         SOR
    };
    
public:
    /*!
     * @param A matrix to build smoother for
     * @param smoother type
     * @param nSweeps number of iterations per smoohting step
     */
    AmrSmoother(const Teuchos::RCP<const matrix_t>& A,
                const Smoother& smoother,
                lo_t nSweeps);
    
    ~AmrSmoother();
    
    /*!
     * Perform one smoothing step
     * @param x right-hand side
     * @param A system matrix
     * @param b right-hand side
     */
    void smooth(const Teuchos::RCP<vector_t>& x,
                const Teuchos::RCP<matrix_t>& A,
                const Teuchos::RCP<vector_t>& b);
    
    /*!
     * Used in AmrMultiGrid constructor
     * @param smoother to create
     */
    static Smoother convertToEnumSmoother(const std::string& smoother);
    
private:
    /*!
     * Initialize paramter list for Ifpack2 smoothers
     * @param smoother of Ifpack2
     * @param nSweeps number of iterations
     */
    void initParameter_m(const Smoother& smoother,
                         lo_t nSweeps);
    
    
private:
    /// Preconditioner instance
    Teuchos::RCP<preconditioner_t> prec_mp;
    
    /// Parameters of preconditioner
    Teuchos::RCP<Teuchos::ParameterList> params_mp;
};

#endif
