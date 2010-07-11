c \BeginDoc
c
c \Name: dnband
c
c \Description:
c
c
c ***************************************************
c * MODIFIED FROM THE ORIGINAL ARPAACK DISTRIBUTION *
c * FOR USE WITH CppNoddy.                          *
c ***************************************************
c This is the original ARPACK dnband wrapper, stripped
c down to only mode 3. This allows it to be called from
c C++ via the cfortran.h interface, with a reduced
c call interface. 
c
c
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style banded form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  dnband can be called with one of the following modes:
c
c  Mode 1:  A*z = lambda*z.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c
c  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*z = amu*z, then 
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
c                 amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*z = lambda*M*z, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*z = amu*z, then 
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c
c  The choice of mode must be specified in IPARAM(7) defined below.
c
c \Usage
c   call dnband
c      ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, 
c        WORKEV, V, N, AB, MB, LDA, RFAC, CFAC, KL, KU, WHICH, 
c        BMAT, NEV, TOL, RESID, NCV, V, LDV, IPARAM, WORKD, 
c        WORKL, LWORKL, WORKC, IWORK, INFO )
c
c \Arguments
c 
c ***  THESE HAVE BEEN TRIMMED DOWN IN THIS REDUCED WRAPPER ***
c
c
c  DR      Double precision array of dimension NEV+1.  (OUTPUT)
c          On exit, DR contains the real part of the Ritz value approximations 
c          to the eigenvalues of A*z = lambda*B*z. 
c
c  DI      Double precision array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR. 
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Real N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit,
c          if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
c          Z represent approximate eigenvectors (Ritz vectors) corresponding
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z computed by DNAUPD.
c
c          The complex Ritz vector associated with the Ritz value
c          with positive imaginary part is stored in two consecutive
c          columns.  The first column holds the real part of the Ritz
c          vector and the second column holds the imaginary part.  The
c          Ritz vector associated with the Ritz value with negative
c          imaginary part is simply the complex conjugate of the Ritz vector
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT) 
c          The leading dimension of the array Z.  If Ritz vectors are 
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.  
c 
c  SIGMAR  Double precision  (INPUT) 
c          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.  
c 
c  SIGMAI  Double precision  (INPUT) 
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the 
c          shift. 
c          Not referenced if IPARAM(7) = 1 or 2.  
c 
c  N       Integer.  (INPUT) 
c          Dimension of the eigenproblem.  
c 
c  AB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC and CFAC. 
c
c  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3 is invoked with a real shift.
c          It is not referenced when IPARAM(7) = 1 or 4.
c
c  CFAC    Complex*16 array of LDA by N. (WORKSPACE/OUTPUT)
c          CFAC is used to store (A-SIGMA*M) and its LU factors
c          when IPARAM(7) = 3 or 4 are used with a complex shift SIGMA.  
c          On exit, it contains the LU factors of (A-SIGMA*M).  
c          It is not referenced when IPARAM(7) = 1 or 2.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c          KU is hardwired to be equal to KL
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c  
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LR' -> want the NEV eigenvalues of largest real part.
c            'SR' -> want the NEV eigenvalues of smallest real part.
c            'LI' -> want the NEV eigenvalues of largest imaginary part.
c            'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c          When IPARAM(7) = 3 or 4, WHICH should be set to 'LM' only. 
c          
c  NEV     Integer. (INPUT)
c          Number of eigenvalues to be computed.
c   
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Arnoldi basis constructed
c          by dnaupd for OP.
c
c  V       Double precision array N by NCV+1.  (OUTPUT)
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c                       represent approximate Schur vectors that span the 
c                       desired invariant subspace.
c          NOTE: The array Z may be set equal to first NEV+1 columns of the 
c          Arnoldi basis vector array V computed by DNAUPD. In this case
c          if RVEC = .TRUE. and HOWMNY='A', then the first NCONV=IPARAM(5) 
c          are the desired Ritz vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: 
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c           ----------------------------------------------------------
c          ISHIFT = 1: exact shift with respect to the current
c                      Hessenberg matrix H.  This is equivalent to
c                      restarting the iteration from the beginning
c                      after updating the starting vector with a linear
c                      combination of Ritz vectors associated with the
c                      "wanted" eigenvalues.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c
c          IPARAM(6) = IUPD
c          Not referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = IPARAM(7): 
c          Must be 3; See under \Description of dnband for the 
c          four modes available.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*z operations,
c                  NUMOPB = total number of B*z operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: The Schur form computed by LAPACK routine dlahqr
c                could not be reordered by LAPACK routine dtrsen.
c                Re-enter subroutine DNEUPD with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: Overflow occurs when we try to transform the Ritz 
c                 values returned from DNAUPD to those of the original
c                 problem using Rayleigh Quotient. 
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
c \EndDoc
c
c------------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
c     May 1995.
c
c\Routines called:
c     dnaupd  ARPACK reverse communication interface routine.
c     dneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine.
c     zgbtrf  LAPACK complex band matrix factorization routine.
c     zgbtrs  LAPACK complex linear system solve routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlamch  LAPACK routine to compute the underflow threshold.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let X' denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the real 
c     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c\Author
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: nband.F   SID: 2.2   DATE OF SID: 11/21/95   RELEASE: 2
c
c\EndLib
c
c---------------------------------------------------------------------
c
c     ORIGINAL CALL INTERFACE
c      subroutine dnband( rvec, howmny, select, dr, di, z, ldz,  sigmar, 
c     &           sigmai, workev, n, ab, mb, lda, rfac,  cfac, kl, ku, 
c     &           which, bmat, nev, tol, resid,  ncv, v, ldv, 
c     &           iparam, workd, workl, lworkl, workc, iwork, info)

c    TRIMMED CALL INTERFACE
      subroutine dnband(dr, di, z, ldz,  sigmar, 
     &           sigmai, n, ab, mb, lda, rfac,  cfac, kl, 
     &           which, nev, resid,  ncv, v, ldv, 
     &           iparam, info)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c 
      character        which*2, bmat, howmny
      integer          n, lda, kl, nev, ncv, ldv,
     &                 ldz, info  
      Double precision sigmar, sigmai 
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*)
      logical          select(ncv)
      Double precision
     &                 dr(*), di(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), rfac(lda,*)
      Complex*16
     &                 cfac(lda,*)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          ipntr(14)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer          ido, i, j, type, imid, itop, ibot, ierr
      Double precision        
     &                 numr, denr, deni, dmdul, safmin 
      logical          rvec, first 
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  one, zero
      parameter        (one = 1.0, zero = 0.0)
c
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
      Double precision
     &                 ddot, dnrm2, dlapy2, dlamch
      external         ddot, dcopy, dgbmv, zgbtrf, zgbtrs, dgbtrf, 
     &                 dgbtrs, dnrm2, dlapy2, dlacpy, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      Intrinsic        dble, dimag, dcmplx
c
c     %-------------------------------------------%
c     | HARDWIRE SOME DETAILS TO TRIM THE WRAPPER |
c     |         REH Sept 16, 2007                 |
c     %-------------------------------------------%
c
      complex*16 workc( n )
      integer iwork( n )
      integer lworkl, ku
      Double precision workev( 3*ncv), workd( 3*n ), 
     &     workl( 3*ncv**2+6*ncv ), tol
      lworkl = 3*ncv**2+6*ncv
      ku = kl
      rvec = .true.
      bmat = 'G'
      howmny = 'A'
      tol = 0.0
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %--------------------------------%
c     | safmin = safe minimum is such  |
c     | that 1/sfmin does not overflow |
c     %--------------------------------%
c
      safmin = dlamch('safmin')
c     
c     %----------------------------------------------------------------%
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     %----------------------------------------------------------------%
c
      type = 4
      if ( iparam(7) .ne. 3 ) then
         print*, ' '
         print*, 'BMAT is inconsistent with IPARAM(7).'
         print*, ' ' 
         go to 9000
      end if
c
c     %------------------------%
c     | Initialize the reverse |
c     | communication flag.    |         
c     %------------------------%
c
      ido   = 0
c
c     %----------------%
c     | Exact shift is |
c     | used.          |
c     %----------------%
c
      iparam(1) = 1
c
c     %-----------------------------------%
c     | Both matrices A and M are stored  |
c     | between rows itop and ibot.  Imid |
c     | is the index of the row that      |
c     | stores the diagonal elements.     |
c     %-----------------------------------%
c
      itop = kl + 1
      imid = kl + ku + 1
      ibot = 2*kl + ku + 1
c
c
c     %-------------------------------------------%
c     | Solving generalized eigenvalue problem in |
c     | shift-invert mode.                        |
c     %-------------------------------------------%
c     
      if ( sigmai .eq. zero ) then
c     
c     %--------------------------------------------%
c     | Construct (A - sigma*M) and factor in real |
c     | arithmetic.                                |
c     %--------------------------------------------%
c     
         do 60 j = 1,n
            do 50 i = itop, ibot 
               rfac(i,j) = ab(i,j) - sigmar*mb(i,j)
 50         continue
 60      continue
c              
         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
         if ( ierr .ne. 0 )  then
            print*, ' '
            print*, '_NBAND: Error with _gbtrf.'
            print*, ' '
            go to 9000
         end if
c     
      else
c     
c     %-----------------------------------------------%
c     | Construct (A - sigma*M) and factor in complex |
c     | arithmetic.                                   |
c     %-----------------------------------------------% 
c     
         do 80 j = 1,n
            do 70 i = itop, ibot 
               cfac(i,j) = dcmplx( ab(i,j)-sigmar*mb(i,j), 
     &              -sigmai*mb(i,j) )
 70         continue 
 80      continue
c     
         call zgbtrf(n, n, kl, ku, cfac, lda, iwork, ierr)
         if ( ierr .NE. 0 )  then
            print*, ' '
            print*, '_NBAND: Error with _gbtrf.'
            print*, ' '
            go to 9000
         end if
c     
      end if 
c end sigmai .eq. zero

c
c     %--------------------------------------------%
c     |  M A I N   L O O P (reverse communication) |
c     %--------------------------------------------%
c
  90  continue 
c
      call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              info )
c
      if (ido .eq. -1) then
c     
c     %-----------------------------------------%
c     | Perform y <-- OP*x                      |
c     |         = Real_part{inv[A-SIGMA*M]*M}*x | 
c     | to force the starting vector into the   |
c     | range of OP.                            |
c     %-----------------------------------------%
c     
         call dgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &        lda, workd(ipntr(1)), 1, zero, 
     &        workd(ipntr(2)), 1)
c     
         if ( sigmai .eq. zero ) then
c     
c     %---------------------%
c     | Shift is real, stay |
c     | in real arithmetic. |
c     %---------------------%            
c     
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &           iwork, workd(ipntr(2)), n, ierr)

            if (ierr .ne. 0) then
               print*, ' ' 
               print*, '_NBAND: Error with _gbtrs.'
               print*, ' ' 
               go to 9000
            end if
c     
         else
c     
c     %--------------------------%
c     | Goto complex arithmetic. |
c     %--------------------------%
c     
            do 120 i = 1,n
               workc(i) = dcmplx(workd(ipntr(2)+i-1))
 120        continue 
c     
            call zgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &           iwork, workc, n, ierr)

            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c     
            do  130 i = 1, n
               workd(ipntr(2)+i-1) = dble(workc(i))
 130        continue 
c     
         end if
c
      else if (ido .eq. 1) then
c     
c     
c     %--------------------------------------%
c     | Perform  y <-- inv(A-sigma*M)*(M*x). |
c     | (M*x) has been computed and stored   |
c     | in workd(ipntr(3)).                  |           
c     %--------------------------------------%
c     
         if ( sigmai .eq. zero ) then
c     
c     %------------------------%
c     | Shift is real, stay in |
c     | real arithmetic.       |
c     %------------------------%
c     
            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)

            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &           iwork, workd(ipntr(2)), n, ierr)

            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c     
         else 
c     
c     %---------------------------%
c     | Go to COMPLEX arithmetic. |
c     %---------------------------%
c     
            do 200 i = 1,n
               workc(i) = dcmplx(workd(ipntr(3)+i-1))
 200        continue 
c     
            call zgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &           iwork, workc, n, ierr)

            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_NBAND: Error in _gbtrs.' 
               print*, ' ' 
               go to 9000
            end if
c     
            do 210 i = 1,n
               workd(ipntr(2)+i-1) = dble(workc(i))
 210        continue 
c     
         end if
c     
      else if (ido .eq. 2) then
c     
c     %--------------------%
c     | Perform y <-- M*x  |
c     %--------------------%
c     
         call dgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &        lda, workd(ipntr(1)), 1, zero, 
     &        workd(ipntr(2)), 1)

c     
      else 
c     
c     %-----------------------------------------%
c     | Either we have convergence, or there is | 
c     | error.                                  |
c     %-----------------------------------------%
c     

         if ( info .lt. 0) then
c     
c     %--------------------------%
c     | Error message, check the |
c     | documentation in DNAUPD  |
c     %--------------------------%
c     
            print *, ' '
            print *, ' Error with _naupd info = ',info
            print *, ' Check the documentation of _naupd '
            print *, ' '
            go to 9000
c     
         else 
c     
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit
     &Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
c     
            if (iparam(5) .gt. 0) then
c     
               call dneupd ( rvec, 'A', select, dr, di, z, ldz, 
     &              sigmar, sigmai, workev, bmat, n, which, 
     &              nev, tol, resid, ncv, v, ldv, iparam,
     &              ipntr, workd, workl, lworkl, info )            
c              

               if ( info .ne. 0) then
c     
c     %------------------------------------%
c     | Check the documentation of DNEUPD. |
c     %------------------------------------%
c     
                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
c     
               else if ( sigmai .ne. zero ) then 
c     
c                  if ( type .eq. 4 ) then
c     
                     first = .true.
                     do 270 j = 1, iparam(5) 
c     
c     %----------------------------------%
c     | Use Rayleigh Quotient to recover |
c     | eigenvalues of the original      |
c     | generalized eigenvalue problem.  |
c     %----------------------------------%
c     
                        if ( di(j) .eq. zero ) then
c     
c     %--------------------------------------%
c     | Eigenvalue is real. Compute          |
c     | d = (x'*inv[A-sigma*M]*M*x) / (x'*x) | 
c     %--------------------------------------%
c     
                           call dgbmv('Nontranspose', n, n, kl, ku, one, 
     $                          mb(itop,1), lda, z(1,j), 1, zero,
     $                          workd, 1)

                           do i = 1, n
                              workc(i) = dcmplx(workd(i))
                           end do
                           call zgbtrs ('Notranspose', n, kl, ku, 1, 
     $                          cfac, lda, iwork, workc, n, info)

                           do i = 1, n
                              workd(i) = dble(workc(i))
                              workd(i+n) = dimag(workc(i))
                           end do
                           denr = ddot(n, z(1,j), 1, workd, 1)
                           deni = ddot(n, z(1,j), 1, workd(n+1), 1)
                           numr  = dnrm2(n, z(1,j), 1)**2
                           dmdul = dlapy2(denr,deni)**2
                           if ( dmdul .ge. safmin ) then
                              dr(j) = sigmar + numr*denr / dmdul
                           else
c     
c     %---------------------%
c     | dmdul is too small. |
c     | Exit to avoid       |
c     | overflow.           |
c     %---------------------%
c     
                              info = -15
                              go to 9000
                           end if
c     
                        else if (first) then 
c     
c     %------------------------%
c     | Eigenvalue is complex. |
c     | Compute the first one  |
c     | of the conjugate pair. |
c     %------------------------%
c     
c     %-------------%
c     | Compute M*x |
c     %-------------%
c     
                           call dgbmv('Nontranspose', n, n, kl, ku,
     $                          one, mb(itop,1), lda, z(1,j), 1, zero,
     $                          workd, 1)

                           call dgbmv('Nontranspose', n, n, kl, ku, 
     $                          one, mb(itop,1), lda, z(1,j+1), 1,
     $                          zero, workd(n+1), 1)

                           do i = 1, n
                              workc(i) = dcmplx(workd(i),workd(i+n))
                           end do
c     
c     %----------------------------%
c     | Compute inv(A-sigma*M)*M*x |
c     %----------------------------%
c     
                           call zgbtrs('Notranspose',n,kl,ku,1,cfac, 
     $                          lda, iwork, workc, n, info)

c     
c     %-------------------------------%
c     | Compute x'*inv(A-sigma*M)*M*x |
c     %-------------------------------%
c     
                           do i = 1, n
                              workd(i) = dble(workc(i))
                              workd(i+n) = dimag(workc(i))
                           end do
                           denr = ddot(n,z(1,j),1,workd,1)
                           denr = denr+ddot(n,z(1,j+1),1,workd(n+1),1)
                           deni = ddot(n,z(1,j),1,workd(n+1),1)
                           deni = deni - ddot(n,z(1,j+1),1,workd,1)
c     
c     %----------------%
c     | Compute (x'*x) |
c     %----------------%
c     
                           numr = dlapy2( dnrm2(n, z(1,j), 1), 
     &                          dnrm2(n, z(1, j+1), 1) )**2
c     
c     %----------------------------------------%
c     | Compute (x'x) / (x'*inv(A-sigma*M)*Mx) |
c     %----------------------------------------%
c     
                           dmdul = dlapy2(denr,deni)**2
                           if ( dmdul .ge. safmin ) then
                              dr(j) = sigmar+numr*denr / dmdul
                              di(j) = sigmai-numr*deni / dmdul
                              first = .false.
                           else
c     
c     %---------------------%
c     | dmdul is too small. |
c     | Exit to avoid       |
c     | overflow.           |
c     %---------------------%
c     
                              info = -15
                              go to 9000
c     
                           end if
c     
                        else
c     
c     %---------------------------%
c     | Get the second eigenvalue |
c     | of the conjugate pair by  |
c     | taking the conjugate of   |
c     | previous one.             |
c     %---------------------------%
c     
                           dr(j) = dr(j-1)
                           di(j) = -di(j-1)
                           first = .true.
c     
                        end if
c     
 270                 continue 
c     
               end if
c     
            end if
c     
         end if
c     
         go to 9000
c     
      end if
c     
c     %----------------------------------------%
c     | L O O P  B A C K to call DNAUPD again. |
c     %----------------------------------------%
c
      go to 90 
c     
 9000 continue
c
      end
