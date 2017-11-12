########## Default compilers & libraries
#
## USER DEFINED COMPONENTS
#
##########

home = "/home/hewitt/CURRENT/Projects/CppNoddy"
#
# Customise these for your compiler/libs/include names/locations
#
c_comp = 'g++'
f_comp = 'gfortran'
#
blas_lib = 'blas'
lapack_lib = 'lapack'
#
# FOR SLEPC/PETSC - makes sure you have the PETSC_DIR, SLEPC_DIR and PETSC_ARCH
# environment variables set. Complex versions of the library must contain the
# string "complex" in $PETSC_ARCH
#
slepc_lib = 'slepc'
petsc_lib = 'petsc'
petsc_inc = '' # we'll construct this using PETSC_DIR and PETSC_ARCH
#
mpi_lib = 'mpi'
mpi_inc = '' # get this for free if you install PETSc with --download-mpich


######################################################
#                                                    #
## YOU PROBABLY SHOULDN'T BE EDITING BELOW HERE      #
#                                                    #
######################################################
#
import os.path
import os

# get the environment
env = Environment(ENV = os.environ)
env['ENV']['TERM'] = os.environ['TERM']

# we need to set the rpath for the linker to find petsc/slepc
# we will populate this below
rpath = []


########## PREAMBLE - get CLI arguments & define colorised output
#
## a slew of command line arguments
#
##########



col = ARGUMENTS.get('col',1)                          # defaults to colourised output
lapack = ARGUMENTS.get('lapack',0)                    # link to LAPACK
slepc = ARGUMENTS.get('slepc',0)	                  # link to SLEPC => PETSC => BLAS/LAPACK
petsc = ARGUMENTS.get('petsc',0)	                  # link to complex PETSC => BLAS/LAPACK (min)
mpi = ARGUMENTS.get('mpi',1)                          # link to MPI
debug = ARGUMENTS.get('debug', 0)                     # ask for debug info to be written to stdout
debug_symbols = ARGUMENTS.get('debug_symbols', 0)     # include debug symbols (-g)
paranoid = ARGUMENTS.get('paranoid', 0)               # paranoid bounds checking
time = ARGUMENTS.get('time', 0)                       # time some routines
profile = ARGUMENTS.get('profile', 0)                 # turn profiling on
private = ARGUMENTS.get('private', 0)                 # compile the private examples?
examples = ARGUMENTS.get('examples', 1)               # compile any of the examples at all?
doc = ARGUMENTS.get('doc', 0)                         # make the documentation

# colour me silly -- unless asked not to with col=0
# <ESC>[{attr};{fg};{bg}m
red = yellow = green = blue = off = ""
if int(col):
    red += "\033[1;31m"
    yellow += "\033[1;33m"
    green += "\033[1;32m"
    blue += "\033[1;34m"
    off += "\033[0m"

# default output format for messaging
def message( col, text ):
    print col + " * " + text + off

#
#
########## DOC BUILD AND EXIT IMMEDIATELY
#
#

# build the documentation if requested then exit
if int(doc):
    print
    message( green, "BUILDING DOCUMENTATION")
    os.system('doxygen')
    message( green, "Exiting")
    Exit(1)

#
#
########## MAIN BUILD STARTS HERE
#
#

try:
    import glob
except ImportError:
    message( red, "No python glob module ... bursting into flames and failing")
    Exit(1)
#
#
########## SETUP OF FILE STRUCTURE
#
#

# source is all cpp files, so let's glob them
src = glob.glob('src/*.cpp')

# set the build dir
topdir = os.getcwd()

# default options
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
libs_str   = 'CppNoddy '
preproc = ' '
opts = ' -O2 -Wall '#-std=c++11 '
link_flags = ' '

#
# extract the PETSc and SLEPc locations from environment variables
petsc_z = petsc_d = False
petsc_dir = os.environ.get('PETSC_DIR','default')
petsc_arch = os.environ.get('PETSC_ARCH','default')
if "complex" in petsc_arch:
    petsc_z = True
    message( red, " PETSC_ARCH points to COMPLEX.")
else:
    petsc_d = True
    message( red, " PETSC_ARCH points to DOUBLE.")
slepc_dir = os.environ.get('SLEPC_DIR','default')

# Set the flags based on command line options
if int(slepc):
    print( " *" )
    message( green, "SLEPc library support is enabled. ")
    message( blue, " This requires PETSc/BLAS/LAPACK support too.")
    message( blue, " Make sure your LD_LIBRARY_PATH contains ")
    if petsc_z:
        message( blue, " $SLEPC_DIR/x86_64-linux-gnu-complex/lib")
        message( blue, " and $PETSC_DIR/x86_64-linux-gnu-complex/lib.")
    else:
        message( blue, " $SLEPC_DIR/x86_64-linux-gnu-double/lib")
        message( blue, " and $PETSC_DIR/x86_64-linux-gnu-double/lib.")
    if (slepc_dir == 'default'):
       message( red, " $SLEPC_DIR is not set!")
       message( red, " You must have $SLEPC_DIR/$PETSC_ARCH point to you SLEPC installation." )
       Exit(1)
    else:
	message( green, " $SLEPC_DIR = "+slepc_dir )
	if ( petsc_arch == 'default' ):
	   message( red, " $PETSC_ARCH is not set!" )
	   message( red, " You must have $PETSC_ARCH be the name of the library build directory in $SLEPC_DIR." )
           Exit(1)
	else:
	   message( green, " $PETSC_ARCH = " + petsc_arch )
	slepc_lib_dir = slepc_dir + "/" + petsc_arch + "/lib"
	slepc_inc = slepc_dir + "/include " + slepc_dir + "/" + petsc_arch + "/include "
	rpath.append( slepc_lib_dir )
    petsc_lib_dir = petsc_dir + "/" + petsc_arch + "/lib"
    petsc_inc = petsc_dir + "/include " + petsc_dir + "/" + petsc_arch + "/include "
    # slepc support requires blas/lapack support
    lapack = 1
    blas = 1
    petsc = 1
    libs_str += slepc_lib + ' ' #+ petsc_lib + ' ' + mpi_lib + ' '
    libdir_str += slepc_lib_dir + ' ' #+ petsc_lib_dir + ' '
    incdir_str += ' ' + slepc_inc + ' ' #+ mpi_inc + ' ' #+ petsc_inc + ' '
    preproc += ' -DSLEPC '

if int(petsc):
    print( " *" )
    message( green, "PETSc library support is enabled. ")
    message( blue, " This requires BLAS/LAPACK/MPI support too.")
    if (petsc_dir == 'default'):
       message( red, " $PETSC_DIR is not set!")
       message( red, " You must have $PETSC_DIR/$PETSC_ARCH point to you PETSC installation." )
       Exit(1)
    else:
        message( green, " $PETSC_DIR = "+petsc_dir )
	if ( petsc_arch == 'default' ):
	   message( red, " $PETSC_ARCH is not set!" )
	   message( red, " You must have $PETSC_ARCH be the name of the library build directory in $SLEPC_DIR." )
	else:
	   message( green, " $PETSC_ARCH = " + petsc_arch )
	petsc_lib_dir = petsc_dir + "/" + petsc_arch + "/lib"
	petsc_inc = petsc_dir + "/include " + petsc_dir + "/" + petsc_arch + "/include "
	rpath.append( petsc_lib_dir )
    # quick hack to determine which PETSc build we have
    if petsc_z:
       preproc += ' -DPETSC_Z '
    else:
       preproc += ' -DPETSC_D '
    # slepc support requires blas/lapack/mumps support
    lapack = 1
    blas = 1
    libs_str += petsc_lib + ' mpi '
    libdir_str += petsc_lib_dir + ' '
    incdir_str += petsc_inc + ' '

if int(lapack):
    print( " *" )
    message( green, "LAPACK support is enabled.")
    message( blue, " I'm assuming a GNU compiler & hard-wiring -DgFortran for")
    message( blue, " cfortran to do its magic with.")
    preproc += ' -DLAPACK -DgFortran'
    libs_str +=  lapack_lib + ' ' + blas_lib + ' '

if int(mpi):
    print( " *" )
    message( green, "MPI library linking is enabled.")
    message( blue, "This is best done via PETSc, and required for MUMPS.")
    libs_str += mpi_lib + ' '
    incdir_str += mpi_inc + ' '
    preproc += ' -DINC_MPI '

if int(debug):
    print( " *" )
    message( red, "DEBUG messages enabled.")
    # overwrite opts to remove any optimisation
    opts = ' '
    debug_symbols = 1
    preproc += ' -DDEBUG '

if int(paranoid):
    print( " *" )
    message( red, "PARANOID checking enabled.")
    preproc += ' -DPARANOID '

if int(time):
    print( " *" )
    message( yellow, "TIMING of selected example problems will be performed.")
    preproc += ' -DTIME '

if int(profile):
    print( " *" )
    message( red, "PROFILING is enabled.")
    opts += ' -pg '
    link_flags += ' -pg '

if int(debug_symbols):
    print( " *" )
    message( yellow, "DEBUG symbols enabled.")
    preproc += ' -g '


print( " *" )
message( blue, " -----  Now checking for dependencies.")

#
#
########## ENVIRONMENT CONSTRUCTION AND CHECKING
#
#


# Split the strings into list elements for lib/inc directoris
incdir = incdir_str.split()
libdir = libdir_str.split()
# we won't split the libs_str just yet, because we might
# change it after checking for libraries


# Construct the environment
env = Environment( FORTRAN = f_comp, CXX = c_comp, CPPPATH = incdir, CCFLAGS = opts + preproc, LINKFLAGS = link_flags, LIBPATH = libdir )

# Now check the environment is complete
conf = Configure( env )

# Check for external libraries if needed
if int(lapack):
    # lapack => blas
    blas = lapack = 1
    # external libraries need fortran support
    if not conf.CheckLib( blas_lib ):
    # check blas
        message( red, "No blas library! It is required to enable LAPACK")
        blas = 0
    if not conf.CheckLib( lapack_lib ):
        # check lapack
        message( red, "No lapack library!")
        lapack = 0
    if ( blas * lapack  == 0 ):
        message( red, "Check your library path ... without the above")
        message( red, "listed libraries, I can't compile with LAPACK support.")
        Exit(1)
    else:
        message( green, "Found BLAS & LAPACK support.")
        message( green, "LAPACK solvers will be used in preference to the native ones.")

if int(petsc):
    if not conf.CheckLib( petsc_lib ):
        message( red, "No libpetsc!")
    else:
	message( green, "Found PETSC.")

if int(slepc):
    if not conf.CheckLib( slepc_lib ):
        message( red, "No libslepc!")
    else:
	message( green, "Found SLEPC.")

if int(slepc):
    # slepc => petsc (and let's assume blas and lapack)
    slepc = blas = lapack = mpi = 1
    if not conf.CheckLib('mpi'):
        message( red, "No libmpi!")
        slepc = 0
    if not conf.CheckLib('slepc'):
        message( red, "No libslepc!")
        slepc = 0
    if not conf.CheckLib('petsc'):
        message( red, "No libpetsc!")
        petsc = 0
    if not conf.CheckLib('lapack'):
        message( red, "No liblapack!")
        lapack = 0
    if not conf.CheckLib('blas'):
        message( red, "No libblas!")
        blas = 0
    if ( slepc * petsc * lapack * blas == 0 ):
        message( red, "SLEPC support has failed.")
        Exit(1)
    else:
        message( green, "Found SLEPc and PETSc, including support for sparse matrix eigensolvers.")


env = conf.Finish()

libs = libs_str.split()

message( blue, " -----  Building.")


#rpath='/home/hewitt/CURRENT/PROGRAMMING/EXTERNAL/petsc-3.7.5/x86_64-linux-gnu-complex/lib'
#
# Build the library in ./lib
#if int(static):
#  print "Static Library"
#  env.StaticLibrary('lib/CppNoddy', src )
#else

print "Shared Library"
env.SharedLibrary('lib/CppNoddy', src, LINKFLAGS = link_flags )

# Call the Examples' SConscript if examples != 0 -- note OPTS & LFLAGS are exported to the build
if int(examples):
    SConscript('Examples/SConscript', exports='env opts preproc link_flags rpath incdir libdir topdir libs private' )
