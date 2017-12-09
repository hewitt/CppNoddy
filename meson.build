project('CppNoddy', 'cpp',
  version : '0.91',
  license : 'GPL2')

#c = run_command('./env.sh')
#message(c.stdout() )

# this is a Cpp project
compiler = meson.get_compiler('cpp')
add_global_arguments( '-O2', language: 'cpp' )

lapack_dep = dependency('lapack', required:false)
blas_dep = dependency('blas', required:false)
petsc_dep = dependency('PETSc', required:false)
slepc_dep = dependency('SLEPc', required:false)
mpich_dep = dependency('mpich', required:false)

slepc = false
petscd = false
petscz = false

if get_option('petscd') or get_option('petscz')
   petsc = true
   if get_option('petscd')
      petscd = true
      message( 'Using real arithmetic PETSc')
   endif
   if get_option('petscz')
      petscz = true
      message( 'Using complex arithmetic PETSc')
   endif
endif
if get_option('slepc')
   slepc = true
   message( 'Using SLEPc')
endif

if lapack_dep.found()
   add_global_arguments('-DLAPACK', language : 'cpp')
endif

if petsc_dep.found()
   if petscd
      add_global_arguments('-DPETSC_D', language : 'cpp')
   endif
   if petscz
      add_global_arguments('-DPETSC_Z', language : 'cpp')
   endif
endif

if slepc_dep.found()
   if slepc
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
                 version : '0.9',
                 name : 'CppNoddy',
                 filebase : 'CppNoddy',
                 description : 'A Noddy numerical computation library.')
