/// \file Exceptions.h
/// The collection of CppNoddy exceptions

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <iostream>
#include <stdexcept>

// Note: instantiation of any of these exceptions will force
// info warning to cout without requiring a what() call if
// WHAT is defined here
#define WHAT

namespace CppNoddy {

  /// An exception to indicate that an error has
  /// been detected in an external (LAPACK) routine
  class ExceptionExternal : public std::runtime_error {
   public:
    ExceptionExternal(const std::string &problem, const int &ifail = 0) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " External library failure \n";
      std::cout << " Returned error code = " << ifail << "\n\n";
#endif

    }

   private:

    void error_header() {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Error: A CppNoddy routine has had an EXTERNAL \n";
      std::cout << "        library problem! \n";
      std::cout << "----------------------------------------------------\n";
    }

  };


  /// An exception class to be thrown when a container
  /// of incorrect geometry used in any class/method.
  class ExceptionGeom : public std::runtime_error {
   public:
    ExceptionGeom(const std::string &problem, const std::size_t &size1, const std::size_t &size2) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " One dimensional container problem \n";
      std::cout << " Container 1 size is " << size1 << "\n";
      std::cout << " Container 2 size is " << size2 << "\n\n";
#endif

    }

    ExceptionGeom(const std::string &problem, const std::size_t &size1, const std::size_t &size2,
                  const std::size_t &size3, const std::size_t &size4) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " Two dimensional container problem \n";
      std::cout << " Container 1 size 1 is " << size1 << "\n";
      std::cout << " Container 1 size 2 is " << size2 << "\n";
      std::cout << " Container 2 size 1 is " << size3 << "\n";
      std::cout << " Container 2 size 2 is " << size4 << "\n\n";
#endif

    }

   private:

    void error_header() {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Error: A CppNoddy routine has had a GEOMETRY error!\n";
      std::cout << "----------------------------------------------------\n";
    }

  };


  /// An exception class that is thrown if too many
  /// Newton steps are taken in either the scalar or
  /// vector Newton classes.
  class ExceptionItn : public std::runtime_error {
   public:
    ExceptionItn(const std::string &problem, const unsigned &itns, const double &max_resid) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " The number of iterations taken so far = " << itns << "\n";
      std::cout << " The maximum residual at this point is = " << max_resid << "\n\n";
#endif

    }

   private:

    void error_header() {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Error: A CppNoddy routine has had an ITERATION error\n";
      std::cout << "----------------------------------------------------\n";
    }

  };


  /// An exception to indicate that a CppNoddy container
  /// has been accessed with index/indices outside the
  /// maximum range for the container.
  class ExceptionRange : public std::runtime_error {
   public:

    ExceptionRange(const std::string &problem, const std::size_t &size1, const std::size_t &index1) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " Range error for one-dimensional container.\n";
      std::cout << " Size 1 of container = " << size1 << "\n";
      std::cout << " Index 1 being accessed = " << index1 << "\n\n";
#endif

    }

    ExceptionRange(const std::string &problem, const std::size_t &size1, const std::size_t &index1,
                   const std::size_t &size2, const std::size_t &index2) :
      std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
      std::cout << " Range error for two-dimensional container.\n";
      std::cout << " Size 1 of container = " << size1 << "\n";
      std::cout << " Size 2 of container = " << size2 << "\n";
      std::cout << " Index 1 being accessed = " << index1 << "\n";
      std::cout << " Index 2 being accessed = " << index2 << "\n\n";
#endif

    }
   private:

    void error_header() {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Error: A CppNoddy routine has had a RANGE error! \n";
      std::cout << "----------------------------------------------------\n";
    }

  };


  /// A generic runtime exception
  class ExceptionRuntime : public std::runtime_error {
   public:
    ExceptionRuntime(const std::string &problem) : std::runtime_error(problem) {
#ifdef WHAT
      error_header();
      std::cout << problem << "\n";
#endif

    }

   private:

    void error_header() {
      std::cout << "----------------------------------------------------\n";
      std::cout << " Error: A CppNoddy routine has had a RUNTIME problem! \n";
      std::cout << "----------------------------------------------------\n";
    }

  };

  /// Not used yet....
  class ExceptionBifurcation {
   public:
    ExceptionBifurcation(const std::string &problem) {
#ifdef WHAT
      std::cout << problem << "\n";
#endif

    }

   private:

  };
}

#endif // EXCEPTIONRUNTIME_H
