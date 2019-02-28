/// \file TwoD_TVDLF_Mesh.h
/// Specification of an object that represents a two dimensional
/// mesh for TVD LF methods.

#ifndef TWOD_TVDLF_MESH_H
#define TWOD_TVDLF_MESH_H

#include <vector>
#include <fstream>
#include <limits>
#include <cassert>

#include <Types.h>
#include <TwoD_TVDLF_Elt.h>
#include <TwoD_Hyperbolic_System.h>

namespace CppNoddy {

  class TwoD_TVDLF_Mesh {
   protected:
    /// iterators for the vector of elements
    typedef std::vector<TwoD_TVDLF_Elt> vector_of_elts;
    typedef vector_of_elts::const_iterator celt_iter;
    typedef vector_of_elts::iterator elt_iter;

    typedef std::vector<TwoD_TVDLF_Elt*> vector_of_boundary_elts;
    typedef vector_of_boundary_elts::iterator bdry_elt_iter;

    /// function pointer used in the initial conditions
    typedef void (*fn_ptr)(const double&, const double&, DenseVector<double>&);

   public:

    /// Constructor for the Finite Volume Mesh using linear elements
    /// \param X A vector of nodal locations at which the element
    ///   FACES will positioned
    /// \param Y A vector of nodal locations at which the element
    ///   FACES will positioned
    /// \param ptr A pointer to the hyperbolic system applied to this mesh
    /// \param init_ptr A pointer to a function that defines the initial conditions
    TwoD_TVDLF_Mesh(const DenseVector<double>& X, const DenseVector<double>& Y,
                    TwoD_Hyperbolic_System* ptr,
                    fn_ptr init_ptr);

    /// Empty desctructor
    virtual ~TwoD_TVDLF_Mesh();

    /// Dump the data to a given filename in a gnuplot format
    /// \param filename The filename to be generated
    void dump_gnu(std::string filename);

    /// Dump the x-nodal positions to a given filename.
    /// This method only applies to regular structured meshes.
    /// \param filename The filename to be generated
    void dump_nodes_x(std::string filename) const;

    /// Dump the y-nodal positions to a given filename.
    /// This method only applies to regular structured meshes.
    /// \param filename The filename to be generated
    void dump_nodes_y(std::string filename) const;

    /// Dump the data over the nodal positions to a given filename.
    /// \param filename The filename to be generated
    void dump_data(std::string filename);

    /// Set the limiter type to be applied in the slope values.
    /// \param id The identifier of the limiter.
    void set_limiter(const unsigned& id);

    /// Update the mesh object. A time step is chosen based upon
    /// the CFL constraint.
    /// \param CFL The CFL constraint for the timestep (e.g. 0.49)
    /// \param max_dt A maximum timestep to be taken irrespective of the CFL value
    /// \return The total timestep that was taken
    double update(const double& CFL, const double& max_dt = std::numeric_limits<long double>::max());

    /// Update the mesh object to a set time level. A time step is chosen based upon
    /// the CFL constraint.
    /// \param CFL The CFL constraint for the timestep (e.g. 0.49)
    /// \param t_end The target end time.
    void update_to(const double& CFL, const double& t_end);

    /// Integrate the concentration values across the entire mesh.
    /// \param mesh_colour The identifier of which mesh to integrate over
    /// \return The vector value of the integral
    DenseVector<double> integrate(std::string mesh_colour = "black");

    /// Given a global coordinate, return a pointer to the elt that contains
    /// that point. Because the elts are stored in a vector, and are logically
    /// in a rectangular grid, we first hunt for the column, then the appropriate
    /// row.
    /// \param x A vector global coordinate
    /// \param mesh_colour A string identifier of the mesh to be searched
    /// \return A pointer to the element containing the global coordinate
    elt_iter get_elt_iter_from_x(const DenseVector<double>& x, std::string mesh_colour = "black");


    /// Get the vector of unknowns at a point in the 2D mesh.
    /// \param x The global coordinate at which to return the values
    /// \return The vector of unknowns
    DenseVector<double> get_point_values(const DenseVector<double>& x);

    /// Get a const reference to the time value for the current mesh.
    /// \return time The time level at which the data in the mesh applies.
    const double& get_time() const;

    /// A virtual method that is run before the first time update.
    /// For user-custom problems.
    /// \param time_step The time step size that is about to be taken
    virtual void actions_before_time_step1(const double& time_step) {
      // Empty by default. If you want to take actions before the time
      // step then you need to inherit from this basic mesh and implement
      // this method.
    }

    /// A virtual method that is run before the second time update.
    /// For user-custom problems.
    /// \param time_step The time step size that is about to be taken
    virtual void actions_before_time_step2(const double& time_step) {
      // Empty by default. If you want to take actions before the time
      // step then you need to inherit from this basic mesh and implement
      // this method.
    }

    /// An STL vector of linear elements -- the black mesh
    vector_of_elts BLACK_ELTS;
    /// An STL vector of linear elements -- the red mesh
    vector_of_elts RED_ELTS;

   protected:

    /// Given an element in the INITIAL STRUCTURED MESH, this
    /// method will return an iterator to an element in the other mesh
    /// that overlaps the corner specified by corner_index.
    /// \param e An element iterator to the source element
    /// \param target_colour The colour of the target mesh, ie. NOT the one that iterator 'e' is in
    /// \param corner_index The index of the corner to be considered
    ///      0,1,2,3 for SW, SE, NE, NW.
    elt_iter get_elt_iter_from_elt_iter(elt_iter e,
                                        std::string target_colour,
                                        int corner_index) {
      // for corners SW,SE,NE,NW:
      // black (i,j) -> red   (i,j),(i+1,j),(i+1,j+1),(i,j+1)
      //  =>  j*Bx+i -> j*Rx+i etc
      // red   (i,j) -> black (i-1,j-1),(i,j-1),(i,j),(i-1,j)
      //  =>  j*Rx+i -> (j-1)*Bx+i-1 etc
      //dimensions of the Black (Bx times By) & Red (Rx times Ry) meshes
      //
      vector_of_elts* target_elts(NULL);
      int k(0);
      int Bx(NX - 1);
      int By(NY - 1);
      int Rx(NX);
      int Ry(NY);
      //
      if(target_colour == "red") {
        int offset(e - BLACK_ELTS.begin());
        target_elts = &RED_ELTS;
        int i, j, l, m;
        //std::cout << " offset =  " << offset << "\n";
        i = offset % Bx;
        j = (offset - i) / Bx;
        switch(corner_index) {
        case 0:
          l = i;
          m = j;
          break;
        case 1:
          l = i + 1;
          m = j;
          break;
        case 2:
          l = i + 1;
          m = j + 1;
          break;
        case 3:
          l = i;
          m = j + 1;
          break;
        default:
          l = 0;
          m = 0;
          assert(false);
        }
        // dont return outside the mesh
        //std::cout << source_colour << " " << i << "," << j << " -> " << l << "," << m <<"\n";
        l = std::min(l, Rx);
        l = std::max(l, 0);
        m = std::min(m, Ry);
        m = std::max(m, 0);
        k = m * Rx + l;
        //std::cout << source_colour << " " << i << "," << j << " -> " << l << "," << m <<"\n";
        //std::cout << "--\n";
      } else if(target_colour == "black") {
        int offset(e - RED_ELTS.begin());
        target_elts = &BLACK_ELTS;
        int i, j, l, m;
        i = offset % Rx;
        j = (offset - i) / Rx;
        switch(corner_index) {
        case 0:
          l = i - 1;
          m = j - 1;
          break;
        case 1:
          l = i;
          m = j - 1;
          break;
        case 2:
          l = i;
          m = j;
          break;
        case 3:
          l = i - 1;
          m = j;
          break;
        default:
          l = 0;
          m = 0;
          assert(false);
        }
        // dont return outside the mesh
        l = std::min(l, Bx);
        l = std::max(l, 0);
        m = std::min(m, By);
        m = std::max(m, 0);
        k = m * Bx + l;
      } else {
        // throw
      }
      return target_elts -> begin() + k;
    }

    /// Given a choice of black or red mesh, return a pointer
    /// to the appropriate vector of elements.
    /// \param mesh_colour A mesh identifier (black or red)
    /// \return A pointer to a std::vector of elements
    vector_of_elts* get_elts_from_colour(std::string mesh_colour);

    /// Given a choice of black or red mesh, return the number
    /// of elements in the x-direction.
    /// \param mesh_colour A mesh identifier (black or red)
    /// \return The number of elts in the x-direction
    std::size_t get_number_elts_in_x(std::string mesh_colour);

    /// Use the appropriate limiter to approximate the slope in each
    /// element in the mesh. The slopes will be sent down to the element
    /// objects and then onto the face objects.
    /// \param elt_vector The vector of elements to set the slope for
    void calc_slopes(vector_of_elts* elt_vector);

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> east_diff(elt_iter e) const;

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> west_diff(elt_iter e) const;

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> north_diff(elt_iter e) const;

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> south_diff(elt_iter e) const;

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> NS_diff(elt_iter e);

    /// Compute a finite difference approximation of the derivative in the
    /// compass direction.
    /// \param e The iterator to the element that we want the derivative for
    /// \return The spatial derivative of the vector system
    DenseVector<double> EW_diff(elt_iter e);

    /// Sign of a double.
    /// \param a The value to return the sign of
    /// \return The sign of the value
    int sgn(double a) const;

    /// A vector version of the minmod operator
    /// \param A A vector to compare
    /// \param B A vector to compare
    /// \return A component-wise minmod vector, i.e., each component
    /// of the returned vector is the minmod of the components of A & B.
    DenseVector<double> minmod(DenseVector<double> A, DenseVector<double> B) const;

    /// A vector version of the maxmod operator
    /// \param A A vector to compare
    /// \param B A vector to compare
    /// \return A component-wise maxmod vector, i.e., each component
    /// of the returned vector is the maxmod of the components of A & B.
    DenseVector<double> maxmod(DenseVector<double> A, DenseVector<double> B) const;

    /// Given an element iterator and the local coordinate this
    /// will return zero for any components specified as inflow
    /// boundary conditions in line with Levy & Tadmor (1997) and
    /// set the centre nodal value to the edge value.
    /// \param e Element iterator
    /// \param face_index The index of the face that is on a boundary
    /// \param diff A vector of derivatives
    void boundary_diff(elt_iter e, const int& face_index, DenseVector<double>& diff) const;

    /// Loops over all boundary elements and sets the Q values in each
    /// one to be the value specified by the edge_values method.
    /// \param elts A vector of elements in the mesh
    void set_boundary_Q(vector_of_elts* elts) {
      // loop over all elts
      elt_iter e(elts -> begin());
      while(e != elts -> end()) {
        // find the external elts that are on the boundary
        if(e -> get_external_flag()) {
          std::set<int> faces(e -> get_external_faces());
          std::set<int>::iterator i(faces.begin());
          // loop through all external faces (there are 2 on corner elts)
          while(i != faces.end()) {
            int face(*i);
            // local coordinate at which to impose the BC
            DenseVector<double> s(2, 0.0);
            switch(face) {
            case 0:
              s[ 0 ] = 0.0;
              s[ 1 ] = -1.0;
              break;
            case 1:
              s[ 0 ] = 1.0;
              s[ 1 ] = 0.0;
              break;
            case 2:
              s[ 0 ] = 0.0;
              s[ 1 ] = 1.0;
              break;
            case 3:
              s[ 0 ] = -1.0;
              s[ 1 ] = 0.0;
              break;
            }
            // get the current edge data
            DenseVector<double> Qe(e -> get_Q(s));
            //// initialise the normal slope vector
            // DenseVector<double> sigma_n( ORDER_OF_SYSTEM, 0.0 );
            // get the user specified edge values
            std::vector<bool> inflow = e -> p_system -> edge_values(face, e -> get_x(s), Qe);     //, sigma_n );
            // get the mid-element nodal data
            DenseVector<double> Qm(e -> get_Q_mid());
            for(std::size_t j = 0; j < ORDER_OF_SYSTEM; ++j) {
              if(inflow[ j ] == true) {
                // only change components that are inflow
                Qm[ j ] = Qe[ j ];
              }
            }
            e -> set_Q_mid(Qm);
            ++i;
          }
        }
        ++e;
      }
    }

    /// Slope limiter method
    unsigned LIMITER;
    /// order of the conservative system
    std::size_t ORDER_OF_SYSTEM;
    /// number of faces in the x & y directions
    std::size_t NX, NY;
    /// the time level of the mesh
    double MESH_TIME;
    /// function pointer to a funnction that defines the initial distribution
    fn_ptr p_Q_INIT;

  };

}

#endif // end OneD_TVDLF_Mesh_H
