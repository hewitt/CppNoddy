project('CppNoddy',
  'cpp',
  version : '0.9',
  license : 'GPL2')

add_global_arguments('-DDEBUG', language : 'cpp')
add_global_arguments('-DPARANOID', language : 'cpp')

lapack_dep = dependency('lapack64','lapack', required:false)
blas_dep = dependency('blas64', 'blas', required:false)
petsc_dep = dependency('petsc', 'Petsc', required:false) # capitalisation matters if not found
slepc_dep = dependency('slepc', 'Slepc', required:false)
mpich_dep = dependency('mpich', required:false)

slepc = false
petscd = false
petscz = false

###############################
# LOOK FOR SPECIFIED OPTIONS  #
###############################

if get_option('petscd')
   petscd = true
   message( 'Using real arithmetic PETSc')
endif
if get_option('petscz')
   petscz = true
   message( 'Using complex arithmetic PETSc')
endif

if get_option('slepc')
   slepc = true
   message( 'Using SLEPc with arithmetic defined by PETSc above')
endif

##########################
# LOOK FOR REQUIRED LIBS #
##########################

if lapack_dep.found()
   message( 'Including LAPACK support' )
   add_global_arguments('-DLAPACK', language : 'cpp')
endif

if petsc_dep.found()
   if petscd
      message( 'Including PETSC_D support' )
      add_global_arguments('-DPETSC_D', language : 'cpp')
   endif
   if petscz
      message( 'Including PETSC_Z support' )
      add_global_arguments('-DPETSC_Z', language : 'cpp')
   endif
endif

if slepc_dep.found()
   if slepc
      message( 'Including SLEPC support' )
      add_global_arguments('-DSLEPC', language : 'cpp')
   endif
endif

# define the CppNoddy includes
inc = include_directories('include')

################################
# SUB DIRECTORIES IN THE BUILD #
################################

subdir('include')
subdir('src')
subdir('Tests')

#######################
# GENERATE THE LIRARY #
#######################

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : CppNoddy,
                 version : '0.9 ',
                 name : 'CppNoddy',
                 filebase : 'CppNoddy',
                 description : 'A Noddy numerical computation library.')
