/// \file SparseLinearEigenSystem.cpp
/// Implementation for the SparseLinearEigenSystem class
/// This class links to SLEPc to perform the solver phase.

#include <vector>
#include <set>
#include <algorithm>
#include <string>

#include <SparseLinearEigenSystem.h>
#include <Exceptions.h>
#include <Types.h>
#include <Timer.h>

#ifdef SLEPC

namespace CppNoddy {

  template <typename _Type>
  SparseLinearEigenSystem<_Type>::SparseLinearEigenSystem(SparseMatrix<_Type >* Aptr, SparseMatrix<_Type >* Bptr) :
    LinearEigenSystem_base() {
    m_region_defined = false;
    m_guess_defined = false;
    //
    m_pA = Aptr;
    m_pB = Bptr;
    m_nev = 8;
    m_nconv = 0;
    m_order = (EPSWhich)7; // target magnitude is the default
    // base class
    m_calc_eigenvectors = true; // SLEPc methods *always* obtain eigenvecs.
  }

  template <typename _Type>
  SparseLinearEigenSystem<_Type>::~SparseLinearEigenSystem() {
  }

  template <typename _Type>
  unsigned SparseLinearEigenSystem<_Type>::get_nconv() const {
    return m_nconv;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::set_nev(unsigned n) {
    m_nev = n;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::set_target(std::complex<double> target) {
    // defaults to (0,0)
    m_shift = target;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::set_order(std::string order_string) {
    int flag(0);
    if(order_string == "EPS_TARGET_MAGNITUDE") {
      m_order=(EPSWhich)7;
      flag=1;
    }
    if(order_string == "EPS_TARGET_REAL") {
      m_order=(EPSWhich)8;
      flag=1;
    }
    if(order_string == "EPS_TARGET_IMAGINARY") {
      m_order=(EPSWhich)9;
      flag=1;
    }
    //typedef enum { EPS_LARGEST_MAGNITUDE=1,
    //EPS_SMALLEST_MAGNITUDE,
    //EPS_LARGEST_REAL,
    //EPS_SMALLEST_REAL,
    //EPS_LARGEST_IMAGINARY,
    //EPS_SMALLEST_IMAGINARY,
    //EPS_TARGET_MAGNITUDE,
    //EPS_TARGET_REAL,
    //EPS_TARGET_IMAGINARY, -- only if COMPLEX
    //EPS_ALL,
    //EPS_WHICH_USER } EPSWhich;
    if(flag==0) {
      std::string problem;
      problem = "The SparseLinearEigenSystem::set_order method has been called\n";
      problem += "with an ordering_string that is not recognised.\n";
      throw ExceptionExternal(problem);
    }
  }

  template <typename _Type>
  bool& SparseLinearEigenSystem<_Type>::region_defined() {
    return m_region_defined;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::set_region(const double& a, const double& b, const double& c, const double& d) {
    m_region_defined = true;
    m_real_left = a;
    m_real_right = b;
    m_imag_bottom = c;
    m_imag_top = d;
  }

  template <typename _Type>
  bool& SparseLinearEigenSystem<_Type>::guess_defined() {
    return m_guess_defined;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::set_initial_guess(const DenseVector<_Type>& guess) {
    m_guess_defined = true;
    m_initial_guess = guess;
  }

  template <typename _Type>
  void SparseLinearEigenSystem<_Type >::eigensolve() {
    // only one method available: SLEPc's generalized shift-invert solver
    eigensolve_slepc();
  }


  template <typename _Type>
  void SparseLinearEigenSystem<_Type>::eigensolve_slepc() {
#ifndef SLEPC
    std::string problem;
    problem = "The SparseLinearEigenSystem::eigensolve method has been called\n";
    problem += "but the compiler option -DSLEPC was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal(problem);
#else
    // create the SLEPc and PETSc data types needed
    Mat petsc_A,petsc_B;
    PetscInt n;
#ifdef PETSC_Z
    // for complex PETSC we only need one (complex) vector to store eigenvec.
    Vec petsc_x;
#endif
#ifdef PETSC_D
    // for double PETSC we need separate real and imag parts for eigenvec.
    Vec petsc_xr,petsc_xi;
#endif
    ST petsc_st;
    EPS petsc_eps;

    // assuming A & B are square
    n = m_pA -> nrows();
    
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Define the matrices that define the eigensystem, Ax=lambdaBx
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // we need to convert from the native sparse format to that required by SLEPc/PETSc

    // define A using PETSc structures
    MatCreate(PETSC_COMM_WORLD,&petsc_A);
    MatSetSizes(petsc_A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(petsc_A);
    
    // we ABSOLUTELY MUST pre-allocate, otherwise the performance really is AWFUL!
    // get the number of non-zero elts in each row as a vector
    PetscInt* all_rows_nnz = new PetscInt[ n ];
    m_pA -> nelts_all_rows(all_rows_nnz);

    // pre-allocate memory using the number of non-zero elts
    // in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(petsc_A, 0, all_rows_nnz);
    // finish the A definition
    MatSetUp(petsc_A);

    
    for(PetscInt i = 0; i<n; ++i) {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = all_rows_nnz[i];
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[all_rows_nnz[i]];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
      // get the data from the CppNoddy sparse matrix structure
      m_pA -> get_row_petsc(i, storage, cols);
      MatSetValues(petsc_A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols;
      delete[] storage;
    }
    // delete the temp storage
    delete[] all_rows_nnz;

    // MatSetValue inserted values are generally cached
    // so we need to explicitly do final assembly
    MatAssemblyBegin(petsc_A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A,MAT_FINAL_ASSEMBLY);

    
    // configure the B matrix
    MatCreate(PETSC_COMM_WORLD,&petsc_B);
    // set B to be an nxn matrix
    MatSetSizes(petsc_B,PETSC_DECIDE,PETSC_DECIDE,n,n);
    // add any command line options
    MatSetFromOptions(petsc_B);

    // we ABSOLUTELY MUST pre-allocate, otherwise the performance really is AWFUL!
    // get the number of non-zero elts in each row as a vector
    all_rows_nnz = new PetscInt[ n ];
    m_pB -> nelts_all_rows(all_rows_nnz);

    // allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(petsc_B, 0, all_rows_nnz);
    // finish the B definition
    MatSetUp(petsc_B);

    // fill the petsc_B matrix from the CppNoddy::SparseMatrix object
    for(PetscInt i = 0; i<n; ++i) {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = all_rows_nnz[i];
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[all_rows_nnz[i]];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
      // get the data from the CppNoddy sparse matrix structure
      m_pB -> get_row_petsc(i, storage, cols);
      MatSetValues(petsc_B,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols;
      delete[] storage;
    }
    // delete the temp storage
    delete[] all_rows_nnz;

    // MatSetValue inserted values are generally cached
    // so we need to explicitly do final assembly
    MatAssemblyBegin(petsc_B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_B,MAT_FINAL_ASSEMBLY);

    // PETSc storage for the eigenvector, using A to define the size
#ifdef PETSC_D
    MatCreateVecs(petsc_A,NULL,&petsc_xr);
    MatCreateVecs(petsc_A,NULL,&petsc_xi);
#endif
#ifdef PETSC_Z
    MatCreateVecs(petsc_A,NULL,&petsc_x);
#endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the eigensolver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // create the eigensolver environment 
    EPSCreate(PETSC_COMM_WORLD,&petsc_eps);
    // define a generalised problem with a B
    EPSSetOperators(petsc_eps,petsc_A,petsc_B);
    // Default is generalizzed non-Hermitian
    EPSSetProblemType(petsc_eps,EPS_GNHEP);
    // Method is Krylov Schur
    EPSSetType(petsc_eps, EPSKRYLOVSCHUR);
    // add any command line options
    EPSSetFromOptions(petsc_eps);
    
    // target spectrum shift - defaults to (0,0)
#ifdef PETSC_D
    EPSSetTarget(petsc_eps, m_shift.real());
#endif
#ifdef PETSC_Z
    EPSSetTarget(petsc_eps, m_shift);
#endif
    // set the order of the returned ev's, as set by the get_order method.
    EPSSetWhichEigenpairs(petsc_eps, m_order);
    // set the number of requested ev's. Not sure if this is relevant if REGION_DEFINED
    //EPSSetDimensions(eps,NEV,NEV+1,PETSC_DEFAULT);
    if ( m_nev < 5 )
    {
      EPSSetDimensions(petsc_eps,m_nev,5,PETSC_DEFAULT);
    }
    else
    {
      EPSSetDimensions(petsc_eps,m_nev,2*m_nev,PETSC_DEFAULT);
    }

    // set tolerance and max number of iterations
    //EPSSetTolerances(petsc_eps, 1.e-6, 20);
    // EPSSetTrueResidual(eps, PETSC_TRUE );
    // EPSSetConvergenceTest(eps, EPS_CONV_ABS);

    // define a monitor function to view convergence (function set above)
    // not needed: use -eps_monitor
    // EPSMonitorSet(petsc_eps,&monitor_function, NULL, NULL);

    //Vec x;
//     VecCreate(PETSC_COMM_WORLD,&petsc_x);
//     VecSetSizes(petsc_x,PETSC_DECIDE,n);
//     VecSetFromOptions(petsc_x);
//     if ( m_guess_defined ) {
//       for ( PetscInt i = 0; i < n; ++i ) {
// #ifdef PETSC_Z
// 	VecSetValue(petsc_x,i,m_initial_guess[i],INSERT_VALUES);
// #endif
// #ifdef PETSC_D
// 	VecSetValue(petsc_x,i,m_initial_guess[i].real(),INSERT_VALUES);
// #endif
//       }
//       #ifdef DEBUG
//       std::cout << "[DEBUG] setting initial space/guess\n";
//       #endif
//       EPSSetInitialSpace(petsc_eps,1,&petsc_x);
//     }
    
    /*
       Define the region containing the eigenvalues of interest
    */
    if(m_region_defined) {
      RG petsc_rg;
      EPSGetRG(petsc_eps, &petsc_rg);
      RGSetType(petsc_rg, RGINTERVAL);
      RGIntervalSetEndpoints(petsc_rg,m_real_left,m_real_right,
                             m_imag_bottom,m_imag_top);
    }

    // get access to the spectral transformation
    EPSGetST(petsc_eps, &petsc_st);
    // we have to use "STSINVERT" instead of the default, because B is
    // typically singular for all problems I'm interested in.
    STSetType(petsc_st, STSINVERT);

    // KSP is the linear solver object of the PETSc library
    KSP petsc_ksp;
    STGetKSP(petsc_st, &petsc_ksp);

    KSPSetOperators(petsc_ksp,petsc_A,petsc_A);
    // set to precondition only
    KSPSetType(petsc_ksp, KSPPREONLY);
    // get a preconditioner object
    PC petsc_pc;
    KSPGetPC(petsc_ksp,&petsc_pc);
    // set it to LU factorization is precondition
    PCSetType(petsc_pc,PCLU);
    PCFactorSetMatSolverType(petsc_pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverType(petsc_pc);

    //  create m_petsc_F
    Mat petsc_F;
    PCFactorGetMatrix(petsc_pc,&petsc_F);

    KSPSetUp(petsc_ksp);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the eigensystem
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    EPSSolve(petsc_eps);
    // EPSGetDimensions(eps,&nev,NULL,NULL);
    // update the NEV private data with the number of returned eigenvalues
    // NEV = (unsigned)nev; // is this always the same as the input nev?

    // Optional: Get some information from the solver and display it
    #ifdef DEBUG
    PetscInt its, lits, maxit;
    EPSType petsc_eps_type;
    PetscReal tol;
    //
    std::cout << "[DEBUG] Target location for eigenvalue  = " << m_shift << "\n";
    std::cout << "[DEBUG] Target ordering of returned eigenvalues (see EPSWhich enum) = " << m_order << "\n";
    #endif
    //
    #ifdef DEBUG
    EPSGetIterationNumber(petsc_eps,&its);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of iterations of the method: %D\n",its);
    EPSGetST(petsc_eps,&petsc_st);
    STGetKSP(petsc_st,&petsc_ksp);
    KSPGetTotalIterations(petsc_ksp,&lits);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of linear iterations of the method: %D\n",lits);
    EPSGetType(petsc_eps,&petsc_eps_type);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Solution method: %s\n\n",petsc_eps_type);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Number of requested eigenvalues: %D\n",m_nev);
    EPSGetTolerances(petsc_eps,&tol,&maxit);
    PetscPrintf(PETSC_COMM_WORLD,"[DEBUG] Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);
    #endif
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
    //EPSReasonView(petsc_eps,PETSC_VIEWER_STDOUT_WORLD); // deprecated in 3.14
    EPSConvergedReasonView(petsc_eps,PETSC_VIEWER_STDOUT_WORLD); // replaces EPSReasonView
    EPSErrorView(petsc_eps,EPS_ERROR_ABSOLUTE,PETSC_VIEWER_STDOUT_WORLD);
    PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
    //#endif

    /*
       Save eigenvectors, if requested
    */
    PetscInt nconv;
    EPSGetConverged(petsc_eps,&nconv);
    // store it in the class
    m_nconv = (unsigned)nconv;
    // create a complex eigenvalue vector
    m_all_eigenvalues = DenseVector<D_complex>(m_nconv, 0.0);
    // complex eigenvector matrix
    m_all_eigenvectors = DenseMatrix<D_complex>(m_nconv, n, 0.0);
    //
    for(unsigned i=0; i<m_nconv; i++) {
#ifdef PETSC_Z
      EPSGetEigenvalue(petsc_eps,i,&m_all_eigenvalues[i],NULL);
      std::cout << m_all_eigenvalues[i] << "\n";
#endif
#ifdef PETSC_D
      double lambda_r,lambda_i;
      EPSGetEigenvalue(petsc_eps,i,&lambda_r,&lambda_i);
      m_all_eigenvalues[i] = D_complex(lambda_r, lambda_i);
#endif

      if(m_calc_eigenvectors) {
#ifdef PETSC_D
        // get the i-th eigenvector from SLEPc
        EPSGetEigenvector(petsc_eps,i,petsc_xr,petsc_xi);
        // convert to a more accessible data structure
        PetscScalar* arrayr;
        VecGetArray1d(petsc_xr, n, 0, &arrayr);
        PetscScalar* arrayi;
        VecGetArray1d(petsc_xi, n, 0, &arrayi);
        for(int j=0; j<n; ++j) {
          m_all_eigenvectors[i][j]=D_complex(arrayr[j], arrayi[j]);
        }
        // documentation says to "restore", though it might not matter as we're done with it now
        VecRestoreArray1d(petsc_xr, n, 0, &arrayr);
        VecRestoreArray1d(petsc_xi, n, 0, &arrayi);
#endif
#ifdef PETSC_Z
        // get the i-th eigenvector from SLEPc
        EPSGetEigenvector(petsc_eps,i,petsc_x,NULL);
        // convert to a more accessible data structure
        PetscScalar* array;
        VecGetArray1d(petsc_x, n, 0, &array);
        for(int j=0; j<n; ++j) {
          m_all_eigenvectors[i][j]=array[j];
        }
        // documentation says to "restore", though it might not matter as we're done with it now
        VecRestoreArray1d(petsc_x, n, 0, &array);
#endif
      }
    }

    /*
       Free work space
    */

    EPSDestroy(&petsc_eps);
    MatDestroy(&petsc_A);
    MatDestroy(&petsc_B);
#ifdef PETSC_D
    VecDestroy(&petsc_xr);
    VecDestroy(&petsc_xi);
#endif
#ifdef PETSC_Z
    VecDestroy(&petsc_x);
#endif

#endif
  }


#ifdef PETSC_Z
  template class SparseLinearEigenSystem<D_complex>
  ;
#endif
#ifdef PETSC_D
  template class SparseLinearEigenSystem<double>
  ;
#endif


} // end namespace

#endif
