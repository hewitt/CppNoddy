/// \file HST.h
/// Some classes useful for hydrodynamic stability theory problems.
/// In particular Cartesian versions of Rayleigh and Orr-Sommerfeld
/// equations.

#ifndef HST_H
#define HST_H

#include <ODE_BVP.h>
#include <Equation_1matrix.h>
#include <DenseLinearEigenSystem.h>

namespace CppNoddy
{
  /// Some utility methods associated with CppNoddy containers.
  namespace HST
  {

    template <typename _Type>
    class Rayleigh
    {
      class Rayleigh_equation : public Equation_1matrix<std::complex<double>, _Type>
      {
      public:
        // this is 3rd order for local refinement
        Rayleigh_equation( ) : Equation_1matrix<std::complex<double> , _Type>( 3 )
        {}

        void residual_fn( const DenseVector<std::complex<double> > &z, DenseVector<std::complex<double> > &g ) const
        {
          _Type y_pos( this -> coord(0) );
          _Type U( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 0 ] );
          _Type Udd( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 1 ] );
          g[ 0 ] = z[ 1 ];
          g[ 1 ] = *p_ALPHA * *p_ALPHA * z[ 0 ] + Udd * z[ 0 ] / ( U - z[ 2 ] );
          g[ 2 ] = 0.0;
        }

        /// The matrix for the BVP coord -- in this case its identity
        void matrix0( const DenseVector<D_complex>& z, DenseMatrix<D_complex>& m ) const
        {
          Utility::fill_identity(m);
        }

        double* p_ALPHA;
        OneD_Node_Mesh<_Type, _Type>* p_BASEFLOW;

      };

      class Rayleigh_left_BC : public Residual<std::complex<double> >
      {
      public:
        // 2 boundary conditions and 3 unknowns
        Rayleigh_left_BC() : Residual<std::complex<double> > ( 2, 3 )
        { }

        void residual_fn( const DenseVector<std::complex<double> > &z, DenseVector<std::complex<double> > &B ) const
        {
          B[ 0 ] = z[ 0 ];
          B[ 1 ] = z[ 1 ] - AMPLITUDE;
        }

        std::complex<double> AMPLITUDE;

      };

      class Rayleigh_right_BC_deriv : public Residual<std::complex<double> >
      {
      public:
        // 1 boundary condition and 3 unknowns
        Rayleigh_right_BC_deriv() : Residual<std::complex<double> > ( 1, 3 ) {}

        void residual_fn( const DenseVector<std::complex<double> > &z, DenseVector<std::complex<double> > &B ) const
        {
          B[ 0 ] = z[ 1 ] + *p_ALPHA * z[ 0 ];
        }

        double* p_ALPHA;
      };

      class Rayleigh_right_BC_Dirichlet : public Residual<std::complex<double> >
      {
      public:
        // 1 boundary condition and 3 unknowns
        Rayleigh_right_BC_Dirichlet() : Residual<std::complex<double> > ( 1, 3 ) {}

        void residual_fn( const DenseVector<std::complex<double> > &z, DenseVector<std::complex<double> > &B ) const
        {
          B[ 0 ] = z[ 0 ];
        }
      };

    public:

      /// ctor -- either for a complex solution in the complex plane,
      /// or a double solution along the real line.
      /// \param base_flow_solution The base flow velocity profile and its curvature.
      /// \param alpha The wavenumber of the perturbation.
      /// \param right_bc_type Determines the choice of boundary condition.
      /// "BL" is a derivative condition, whilst "CHANNEL" is a Dirichlet
      /// impermeability condition

      Rayleigh( OneD_Node_Mesh<_Type, _Type> &base_flow_solution, double& alpha, const std::string &right_bc_type = "BL" );

      /// Solve the global eigenvalue problem for the Rayleigh equation
      /// by employing a second-order finite-difference matrix. The
      /// scheme allows for a non-uniform distribution of nodal points.
      void global_evp( );

      /// Solve the EVP locally as a nonlinear BVP for just one mode.
      /// Again we employ a second-order accurate finite-difference scheme
      /// which allows for non-uniform distribution of nodal points.
      /// \param i_ev The index of the eigenvalue to solve for, based on the return
      /// from the global_evp method.
      void local_evp( std::size_t i_ev );

      /// Refine the EIGENVECTORS mesh for a new baseflow. Useful following
      /// a global_evp solve and prior to a local_evp solve.
      /// \param new_baseflow A new mesh containing the base flow, we'll linearly
      /// interpolate to this mesh
      void remesh1( const OneD_Node_Mesh<_Type, _Type>& new_baseflow );

      /// A handle to the eigenvectors mesh
      /// \return A mesh containing the eigenvectors, with the i-th variable
      /// corresponding to the i-th eigenvalue.
      OneD_Node_Mesh<std::complex<double> , _Type>& eigenvectors()
      {
        return EIGENVECTORS;
      }

      OneD_Node_Mesh<std::complex<double> > eigenvector( unsigned i_ev )
      {
        unsigned n( BASEFLOW.get_nnodes() );
        OneD_Node_Mesh<std::complex<double> > temp( BASEFLOW.nodes(), 1 );
        for ( unsigned i = 0; i < n; ++i )
        {
          temp( i, 0 ) = EIGENVECTORS( i, i_ev );
        }
        return temp;
      }

      /// Iterate on the wavenumber ALPHA, using the local_evp routine, to drive a
      /// selected eigenvalue to be neutral (ie. imaginary part is zero)
      /// \param i_ev The index of the eigenvalue to iterate on
      void iterate_to_neutral( std::size_t i_ev );

      /// A handle to the eigenvalues vector
      /// \return A vector of eigenvalues
      DenseVector<std::complex<double> >& eigenvalues()
      {
        return EIGENVALUES;
      }

      /// A handle to the wavenumber
      /// \return The wavenumber
      double& alpha()
      {
        return ALPHA;
      }

    protected:

      std::string RIGHT_BC_TYPE;
      double ALPHA;
      OneD_Node_Mesh<_Type, _Type> BASEFLOW;
      OneD_Node_Mesh<std::complex<double> , _Type> EIGENVECTORS;
      DenseVector<std::complex<double> > EIGENVALUES;
    };




    class Orr_Sommerfeld
    {
      enum { phi, phid, psi, psid, eval };

      class Orr_Sommerfeld_equation : public Equation_1matrix<std::complex<double> >
      {
      public:
        // this is 5th order for local refinement
        Orr_Sommerfeld_equation( ) : Equation_1matrix<std::complex<double> >( 5 )
        {}

        void residual_fn( const DenseVector<std::complex<double> > &z, DenseVector<std::complex<double> > &g ) const
        {
          double y_pos( this -> coord(0) );
          double U( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 0 ] );
          double Udd( p_BASEFLOW -> get_interpolated_vars( y_pos )[ 1 ] );
          // define the equation as 5 1st order equations
          g[ phi ] = z[ phid ];
          g[ phid ] = z[ psi ] + *p_ALPHA * *p_ALPHA * z[ phi ];
          g[ psi ] = z[ psid ];
          g[ psid ] = *p_ALPHA * *p_ALPHA * z[ psi ]
                      + D_complex( 0.0, 1.0 ) * *p_ALPHA * *p_RE * ( U * z[ psi ] - Udd * z[ phi ] )
                      - D_complex( 0.0, 1.0 ) * *p_ALPHA * *p_RE * z[ eval ] * z[ psi ];
          g[ eval ] = 0.0;
        }

        /// The matrix for the BVP coord -- in this case it's an identity matrix
        void matrix0( const DenseVector<D_complex>& z, DenseMatrix<D_complex>& m ) const
        {
          Utility::fill_identity(m);
        }

        double* p_ALPHA;
        double* p_RE;
        OneD_Node_Mesh<double>* p_BASEFLOW;

      };

      class OS_left_BC : public Residual<D_complex>
      {
      public:
        // 3 boundary conditions and 5 unknowns
        OS_left_BC() : Residual<D_complex> ( 3, 5 ) {}

        void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
        {
          B[ 0 ] = z[ phi ];
          B[ 1 ] = z[ phid ];
          B[ 2 ] = z[ psi ] - 1.0; // an arbitrary amplitude traded against the eigenvalue
        }
      };

      class OS_right_BC : public Residual<D_complex>
      {
      public:
        // 2 boundary conditions and 5 unknowns
        OS_right_BC() : Residual<D_complex> ( 2, 5 ) {}

        void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
        {
          B[ 0 ] = z[ phi ];
          B[ 1 ] = z[ phid ];
        }
      };


    public:

      /// ctor
      /// \param base_flow_solution The base flow velocity profile and
      /// its curvature.
      /// \param alpha The wavenumber to compute the spectrum for.
      /// \param rey The Reynolds number to compute the spectrum for.
      Orr_Sommerfeld( OneD_Node_Mesh<double> &base_flow_solution, double alpha, double rey )
      {
        ALPHA = alpha;
        RE = rey;
        BASEFLOW = base_flow_solution;
      }

      /// Solve the global eigenvalue problem for the Rayleigh equation.
      void global_evp( )
      {
        const std::complex<double> I( 0.0, 1.0 );
        unsigned N( 2 * BASEFLOW.get_nnodes() );
        // matrices for the EVP, initialised with zeroes
        DenseMatrix<D_complex> a( N, N, 0.0 );
        DenseMatrix<D_complex> b( N, N, 0.0 );
        //
        double d = BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 );
        // boundary conditions at the left boundary
        a( 0, 0 ) = 1.0;           // phi( left ) = 0
        a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
        a( 1, 2 ) = 2.0 / d;
        a( 1, 4 ) = -0.5 / d;
        // fill the interior nodes
        for ( std::size_t i = 1; i <= BASEFLOW.get_nnodes() - 2; ++i )
        {
          // position in the channel
          //const double y = BASEFLOW.coord( i );
          // base flow profile
          const double U = BASEFLOW(i,0);
          // BASEFLOW.get_interpolated_vars( y )[ 0 ];
          const double Udd = (BASEFLOW(i+1,0) - 2*BASEFLOW(i,0) + BASEFLOW(i-1,0))/(d*d);
          //BASEFLOW.get_interpolated_vars( y )[ 1 ];

          // the first quation at the i'th nodal point
          std::size_t row = 2 * i;
          a( row, row ) = -2.0 / ( d * d ) - ALPHA * ALPHA;
          a( row, row - 2 ) = 1.0 / ( d * d );
          a( row, row + 2 ) = 1.0 / ( d * d );
          a( row, row + 1 ) = -1.0;

          row += 1;
          // the second equation at the i'th nodal point
          a( row, row ) = -2.0 / ( d * d ) - ALPHA * ALPHA - I * ALPHA * RE * U;
          a( row, row - 2 ) = 1.0 / ( d * d );
          a( row, row + 2 ) = 1.0 / ( d * d );
          a( row, row - 1 ) = I * ALPHA * RE * Udd;

          b( row, row ) = - I * ALPHA * RE;
        }
        // boundary conditions at right boundary
        a( N - 2, N - 2 ) = 1.5 / d;
        a( N - 2, N - 4 ) = -2.0 / d;
        a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
        a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0
        // a vector for the eigenvalues
        DenseLinearEigenSystem<D_complex> system( &a, &b );
        system.eigensolve();
        system.tag_eigenvalues_disc( +1, 100.0 );
        EIGENVALUES = system.get_tagged_eigenvalues();

      }

      /// Solve the EVP locally as a nonlinear BVP for just one mode.
      /// \param i_ev The index of the eigenvalue to solve for, based on the return
      /// from the global_evp method.
      void local_evp( std::size_t i_ev ) {}

      /// Refine the EIGENVECTORS mesh for a new baseflow. Useful following
      /// a global_evp solve and prior to a local_evp solve.
      /// \param new_baseflow A new mesh containing the base flow, we'll linearly
      /// interpolate to this mesh
      void remesh1( const OneD_Node_Mesh<double>& new_baseflow ) {}

      /// A handle to the eigenvectors mesh
      /// \return A mesh containing the eigenvectors, with the i-th variable
      /// corresponding to the i-th eigenvalue.
      OneD_Node_Mesh<std::complex<double> >& eigenvectors()
      {
        return EIGENVECTORS;
      }


      /// Iterate on the wavenumber ALPHA, using the local_evp routine, to drive a
      /// selected eigenvalue to be neutral (ie. imaginary part is zero)
      /// \param i_ev The index of the eigenvalue to iterate on
      void iterate_to_neutral( std::size_t i_ev ) {}

      /// A handle to the eigenvalues vector
      /// \return A vector of eigenvalues
      DenseVector<std::complex<double> >& eigenvalues()
      {
        return EIGENVALUES;
      }

      /// A handle to the wavenumber
      /// \return The wavenumber
      double& alpha()
      {
        return ALPHA;
      }

    protected:

      double RE;
      double ALPHA;
      OneD_Node_Mesh<double > BASEFLOW;
      OneD_Node_Mesh<std::complex<double> > EIGENVECTORS;
      DenseVector<std::complex<double> > EIGENVALUES;
    };







    template <typename _Type>
    Rayleigh<_Type>::Rayleigh( OneD_Node_Mesh<_Type, _Type> &base_flow_solution, double& alpha, const std::string &right_bc_type )
    {
      // protected data storage
      RIGHT_BC_TYPE = right_bc_type;
      ALPHA = alpha;
      BASEFLOW = base_flow_solution;
      EIGENVALUES = DenseVector<std::complex<double> >( BASEFLOW.get_nnodes(), 0.0 );
      EIGENVECTORS = OneD_Node_Mesh<std::complex<double>, _Type>( BASEFLOW.nodes(), 1 );
    }


    template <typename _Type>
    void Rayleigh<_Type>::remesh1( const OneD_Node_Mesh<_Type, _Type>& new_baseflow )
    {
      // reset the baseflow mesh
      BASEFLOW = new_baseflow;
      // first order interpolation of the eigenvectors
      EIGENVECTORS.remesh1( BASEFLOW.nodes() );
    }

    template<>
    void Rayleigh<double>::iterate_to_neutral( std::size_t i_ev ) {}

    template<>
    void Rayleigh<std::complex<double> >::iterate_to_neutral( std::size_t i_ev )
    {
      double delta( 1.e-8 );
      do
      {
        std::cout << "ALPHA = " << ALPHA << "\n";
        local_evp( i_ev );
        std::complex<double> copy_of_ev( EIGENVALUES[ i_ev ] );
        ALPHA += delta;
        local_evp( i_ev );
        ALPHA -= delta;
        double d_ev = ( std::imag( EIGENVALUES[ i_ev ] ) - std::imag( copy_of_ev ) ) / delta;
        ALPHA -= std::imag( copy_of_ev ) / d_ev;
        std::cout << "ITERATING: " << ALPHA << " " << EIGENVALUES[ i_ev ] << " " << d_ev << "\n";
      }
      while ( std::abs( std::imag( EIGENVALUES[ i_ev ] ) ) > 1.e-6 );
    }


    template <typename _Type>
    void Rayleigh<_Type>::local_evp( std::size_t i_ev  )
    {
      // number of nodes in the mesh
      std::size_t N = BASEFLOW.get_nnodes();
      // formulate the Rayleigh equation as a BVP
      Rayleigh_equation Rayleigh_problem;
      // boundary conditions
      Rayleigh_left_BC Rayleigh_left;
      Rayleigh_right_BC_deriv Rayleigh_right_deriv;
      Rayleigh_right_BC_Dirichlet Rayleigh_right_Dirichlet;
      // set the private member data in the objects
      Rayleigh_problem.p_BASEFLOW = &BASEFLOW;
      Rayleigh_problem.p_ALPHA = &ALPHA;
      Rayleigh_right_deriv.p_ALPHA = &ALPHA;

      // pointer to the equation
      ODE_BVP<std::complex<double>, _Type>* p_Rayleigh;
      if ( RIGHT_BC_TYPE == "BL" )
      {
        p_Rayleigh = new ODE_BVP<std::complex<double>, _Type>( &Rayleigh_problem, BASEFLOW.nodes(), &Rayleigh_left, &Rayleigh_right_deriv );
      }
      else
        if ( RIGHT_BC_TYPE == "CHANNEL" )
        {
          p_Rayleigh = new ODE_BVP<std::complex<double>, _Type>( &Rayleigh_problem, BASEFLOW.nodes(), &Rayleigh_left, &Rayleigh_right_Dirichlet );
        }
        else
        {
          std::string problem;
          problem = " The Rayleigh.global_evp_uniform class has been called with an unknown string \n";
          problem += " that defines the far-field boundary condition.";
          throw ExceptionRuntime( problem );
        }

      p_Rayleigh -> max_itns() = 30;

      // set the initial guess using the global_evp solve data
      p_Rayleigh -> solution()( 0, 0 ) = EIGENVECTORS( 0, i_ev );
      p_Rayleigh -> solution()( 0, 1 ) = ( EIGENVECTORS( 1, i_ev ) - EIGENVECTORS( 0, i_ev ) ) / ( BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 ) );
      p_Rayleigh -> solution()( 0, 2 ) = EIGENVALUES[ i_ev ];
      // set the (arbitrary) amplitude from the global_evp solve data
      Rayleigh_left.AMPLITUDE = p_Rayleigh -> solution()( 0, 1 );
      for ( unsigned i = 1; i < N - 1; ++i )
      {
        p_Rayleigh -> solution()( i, 0 ) = EIGENVECTORS( i, i_ev );
        p_Rayleigh -> solution()( i, 1 ) = ( EIGENVECTORS( i + 1, i_ev ) - EIGENVECTORS( i - 1, i_ev ) ) / ( BASEFLOW.coord( i + 1 ) - BASEFLOW.coord( i - 1 ) );
        p_Rayleigh -> solution()( i, 2 ) = EIGENVALUES[ i_ev ];
      }
      p_Rayleigh -> solution()( N - 1, 0 ) = EIGENVECTORS( N - 1, i_ev );
      p_Rayleigh -> solution()( N - 1, 1 ) = ( EIGENVECTORS( N - 1, i_ev ) - EIGENVECTORS( N - 2, i_ev ) ) / ( BASEFLOW.coord( N - 1 ) - BASEFLOW.coord( N - 2 ) );
      p_Rayleigh -> solution()( N - 1, 2 ) = EIGENVALUES[ i_ev ];

      // do a local solve
      p_Rayleigh -> solve2();
      // write the eigenvalue and eigenvector back to private member data store
      EIGENVALUES[ i_ev ] = p_Rayleigh -> solution()( 0, 2 );
      for ( unsigned i = 0; i < N; ++i )
      {
        EIGENVECTORS( i, i_ev ) = p_Rayleigh -> solution()( i, 0 );
      }
      // delete the equation object
      delete p_Rayleigh;
    }


    template <>
    void Rayleigh<double>::global_evp( )
    {
      unsigned N( BASEFLOW.get_nnodes() );
      // Rayleigh EVP -- we'll keep it complex even though its extra work
      DenseMatrix<std::complex<double> > A( N, N, 0.0 );
      DenseMatrix<std::complex<double> > B( N, N, 0.0 );
      // f_0 = 0
      A( 0, 0 ) = 1.0;
      B( 0, 0 ) = 0.0;
      // step through the interior nodes
      for ( std::size_t i = 1; i < N - 1; ++i )
      {
#ifdef PARANOID
        const double h1 = BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 );
        if ( std::abs( BASEFLOW.coord( i ) - BASEFLOW.coord( i - 1 ) - h1 ) > 1.e-12 )
        {
          std::string problem;
          problem = " Warning: The Rayleigh.global_evp method is only first order \n";
          problem += " on a non-uniform mesh.\n";
          throw ExceptionRuntime( problem );
        }
#endif
        double h( BASEFLOW.coord( i ) - BASEFLOW.coord( i - 1 ) );
        double k( BASEFLOW.coord( i + 1 ) - BASEFLOW.coord( i ) );
        double sigma( k / h ); // sigma = 1 => uniform mesh
        double y = BASEFLOW.coord( i );
        double U = BASEFLOW.get_interpolated_vars( y )[ 0 ];
        double Udd = BASEFLOW.get_interpolated_vars( y )[ 1 ];
        double h2 = 0.5 * h * h * sigma * ( sigma + 1.0 );
        A( i, i ) = U * ( -( sigma + 1.0 ) / h2 - ALPHA * ALPHA ) - Udd;
        A( i, i - 1 ) = sigma * U / h2;
        A( i, i + 1 ) = U / h2;
        //
        B( i, i ) = - ( sigma + 1.0 ) / h2 - ALPHA * ALPHA;
        B( i, i - 1 ) = sigma / h2;
        B( i, i + 1 ) = 1. / h2;
      }

      // 3 point backward difference the far-field in BL, but pin in channel
      if ( RIGHT_BC_TYPE == "BL" )
      {
        double h = BASEFLOW.coord( N - 2 ) - BASEFLOW.coord( N - 3 );
        double k = BASEFLOW.coord( N - 1 ) - BASEFLOW.coord( N - 2 );
        //A( N - 1, N - 1 ) = -3.0 / ( 2 * h )- alpha;
        //A( N - 1, N - 2 ) = 4.0 / ( 2 * h );
        //A( N - 1, N - 3 ) = -1.0 / ( 2 * h );
        A( N - 1, N - 1 ) = h * ( h + k ) / k * ( - 1. / ( h * h ) + 1. / ( ( h + k ) * ( h + k ) ) ) - ALPHA;
        A( N - 1, N - 2 ) = h * ( h + k ) / k * ( 1. / ( h * h ) );
        A( N - 1, N - 3 ) = h * ( h + k ) / k * ( 1. / ( ( h + k ) * ( h + k ) ) );
        B( N - 1, N - 1 ) = 0.0;
      }
      else
        if ( RIGHT_BC_TYPE == "CHANNEL" )
        {
          A( N - 1, N - 1 ) = 1.0;
          B( N - 1, N - 1 ) = 0.0;
        }
        else
        {
          std::string problem;
          problem = " The Rayleigh.global_evp_uniform class has been called with an unknown string \n";
          problem += " that defines the far-field boundary condition.";
          throw ExceptionRuntime( problem );
        }
      double U_max( BASEFLOW( 0, 0 ) );
      double U_min( BASEFLOW( 0, 0 ) );
      for ( unsigned i = 1; i < N; ++i )
      {
        U_max = std::max( U_max, BASEFLOW( i, 0 ) );
        U_min = std::min( U_min, BASEFLOW( i, 0 ) );
      }

      DenseLinearEigenSystem<std::complex<double> > rayleigh_evp( &A, &B );
      rayleigh_evp.eigensolve();
      // return based on Howard's circle theorem
      rayleigh_evp.set_shift( std::complex<double> ( 0.5 * ( U_max + U_min ), 0.0 ) );
      rayleigh_evp.tag_eigenvalues_disc( +1, 0.5 * ( U_max - U_min ) );
      EIGENVALUES = rayleigh_evp.get_tagged_eigenvalues();
      DenseMatrix<std::complex<double> > eigenvecs_mtx = rayleigh_evp.get_tagged_eigenvectors();
      // get eigenvectors
      EIGENVECTORS = OneD_Node_Mesh<std::complex<double>, double>( BASEFLOW.nodes(), N );
      for ( unsigned evec = 0; evec < EIGENVALUES.size(); ++evec )
      {
        for ( unsigned node = 0; node < N; ++node )
        {
          EIGENVECTORS( node, evec ) = eigenvecs_mtx( evec, node );
        }
      }
    }


    template <>
    void Rayleigh<std::complex<double> >::global_evp( )
    {
      unsigned N( BASEFLOW.get_nnodes() );
      // Rayleigh EVP
      DenseMatrix<std::complex<double> > A( N, N, 0.0 );
      DenseMatrix<std::complex<double> > B( N, N, 0.0 );
      // f_0 = 0
      A( 0, 0 ) = 1.0;
      B( 0, 0 ) = 0.0;
      //const std::complex<double>  h1 = BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 );
      // step through the interior nodes
      for ( std::size_t i = 1; i < N - 1; ++i )
      {
#ifdef PARANOID
        const double h1 = abs(BASEFLOW.coord( 1 ) - BASEFLOW.coord( 0 ));
        if ( std::abs( BASEFLOW.coord( i ) - BASEFLOW.coord( i - 1 ) - h1 ) > 1.e-12 )
        {
          std::string problem;
          problem = " Warning: The Rayleigh.global_evp method is only first order \n";
          problem += " on a non-uniform mesh.\n";
          //throw ExceptionRuntime( problem );
        }
#endif
        std::complex<double>  h( BASEFLOW.coord( i ) - BASEFLOW.coord( i - 1 ) );
        std::complex<double>  k( BASEFLOW.coord( i + 1 ) - BASEFLOW.coord( i ) );
        std::complex<double>  sigma( k / h ); // sigma = 1 => uniform mesh
        std::complex<double>  y = BASEFLOW.coord( i );
        std::complex<double>  U = BASEFLOW.get_interpolated_vars( y )[ 0 ];
        std::complex<double>  Udd = BASEFLOW.get_interpolated_vars( y )[ 1 ];
        std::complex<double>  h2 = 0.5 * h * h * sigma * ( sigma + 1.0 );
        A( i, i ) = U * ( -( sigma + 1.0 ) / h2 - ALPHA * ALPHA ) - Udd;
        A( i, i - 1 ) = sigma * U / h2;
        A( i, i + 1 ) = U / h2;
        //
        B( i, i ) = - ( sigma + 1.0 ) / h2 - ALPHA * ALPHA;
        B( i, i - 1 ) = sigma / h2;
        B( i, i + 1 ) = 1. / h2;
      }

      // 3 point backward difference the far-field in BL, but pin in channel
      if ( RIGHT_BC_TYPE == "BL" )
      {
        std::complex<double>  h = BASEFLOW.coord( N - 2 ) - BASEFLOW.coord( N - 3 );
        std::complex<double>  k = BASEFLOW.coord( N - 1 ) - BASEFLOW.coord( N - 2 );
        //A( N - 1, N - 1 ) = -3.0 / ( 2 * h )- alpha;
        //A( N - 1, N - 2 ) = 4.0 / ( 2 * h );
        //A( N - 1, N - 3 ) = -1.0 / ( 2 * h );
        A( N - 1, N - 1 ) = h * ( h + k ) / k * ( - 1. / ( h * h ) + 1. / ( ( h + k ) * ( h + k ) ) ) - ALPHA;
        A( N - 1, N - 2 ) = h * ( h + k ) / k * ( 1. / ( h * h ) );
        A( N - 1, N - 3 ) = h * ( h + k ) / k * ( 1. / ( ( h + k ) * ( h + k ) ) );
        B( N - 1, N - 1 ) = 0.0;
      }
      else
        if ( RIGHT_BC_TYPE == "CHANNEL" )
        {
          A( N - 1, N - 1 ) = 1.0;
          B( N - 1, N - 1 ) = 0.0;
        }
        else
        {
          std::string problem;
          problem = " The Rayleigh.global_evp_uniform class has been called with an unknown string \n";
          problem += " that defines the far-field boundary condition.";
          throw ExceptionRuntime( problem );
        }
      double U_max( BASEFLOW( 0, 0 ).real() );
      double U_min( BASEFLOW( 0, 0 ).real() );
      for ( unsigned i = 1; i < N; ++i )
      {
        U_max = std::max( U_max, BASEFLOW( i, 0 ).real() );
        U_min = std::min( U_min, BASEFLOW( i, 0 ).real() );
      }

      DenseLinearEigenSystem<std::complex<double> > rayleigh_evp( &A, &B );
      rayleigh_evp.eigensolve();
      // return based on Howard's circle theorem
      rayleigh_evp.set_shift( std::complex<double> ( 0.5 * ( U_max + U_min ), 0.0 ) );
      rayleigh_evp.tag_eigenvalues_disc( +1, 0.5 * ( U_max - U_min ) );
      EIGENVALUES = rayleigh_evp.get_tagged_eigenvalues();
      DenseMatrix<std::complex<double> > eigenvecs_mtx = rayleigh_evp.get_tagged_eigenvectors();
      // get eigenvectors
      EIGENVECTORS = OneD_Node_Mesh<std::complex<double>, std::complex<double> >( BASEFLOW.nodes(), N );
      for ( unsigned evec = 0; evec < EIGENVALUES.size(); ++evec )
      {
        for ( unsigned node = 0; node < N; ++node )
        {
          EIGENVECTORS( node, evec ) = eigenvecs_mtx( evec, node );
        }
      }
    }

  }
} // end namespace

#endif // HST_H
