########## Default compilers & libraries
#
## USER DEFINED COMPONENTS
#
##########

# Customise these for your compiler/libs/include names/locations
c_comp = 'g++'
f_comp = 'gfortran'
blas_lib = 'blas'
lapack_lib = 'lapack'
superlu_lib = 'superlu'
superlu_inc = '/usr/include/superlu'

#
#
########## YOU PROBABLY SHOULDN'T BE EDITING BELOW HERE
#
#

import os.path

########## PREAMBLE - get CLI arguments & define colorised output
#
## a slew of command line arguments
#
##########

col = ARGUMENTS.get('col',1)                          # defaults to colourised output
static = ARGUMENTS.get('static',0)                    # static library
lapack = ARGUMENTS.get('lapack',0)                    # link to LAPACK
arpack = ARGUMENTS.get('arpack',0)                    # link to ARPACK -- just a test case, do not use
superlu = ARGUMENTS.get('superlu',0)                  # link to SUPERLU
debug = ARGUMENTS.get('debug', 0)                     # ask for debug info to be written to stdout
debug_symbols = ARGUMENTS.get('debug_symbols', 0)     # include debug symbols (-g)
paranoid = ARGUMENTS.get('paranoid', 0)               # paranoid bounds checking
warn = ARGUMENTS.get('warn', 1)                       # compilation warnings on (defaults to on)
time = ARGUMENTS.get('time', 0)                       # time some routines
profile = ARGUMENTS.get('profile', 0)                 # turn profiling on
private = ARGUMENTS.get('private', 0)                 # compile the private examples?
examples = ARGUMENTS.get('examples', 1)               # compile any of the examples at all?
doc = ARGUMENTS.get('doc', 0)                         # make the documentation
force_g2c = ARGUMENTS.get('force_g2c', 0)             # force use of g2c instead of gfortran (legacy)

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
    message( red, "No globbing module ... bursting into flames and failing")
    Exit(1)
#
#
########## SETUP OF FILE STRUCTURE
#
#

# source is all cpp files, so let's glob them
# some fortran will be added later for arpack wrapping-up
src = glob.glob('src/*.cpp')

# set the build dir
topdir = os.getcwd()

# default options
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
libs_str   = 'CppNoddy '
preproc = ' '
opts = '-O3 ' 
link_flags = ' '


# Set the flags based on command line options
if int(arpack):
    message( green, "ARPACK support is enabled. This require BLAS/LAPACK support too.")
    message( red, "I'm assuming a GNU compiler & hard-wiring -DgFortran for")
    message( red, "cfortran to do its magic with.")
    # arpack support requires lapack support
    lapack = 1
    libs_str += 'arpack gfortran '
    preproc += ' -DARPACK'
    # add the fortran banded arpack wrappers to the source suite
    src +=  glob.glob('src/*.f')

if int(superlu):
    message( green, "SUPERLU support is enabled. This requires BLAS support too.")
    # arpack support requires blas support
    lapack = 1
    libs_str += superlu_lib + ' ' + blas_lib + ' '
    incdir_str += ' ' + superlu_inc + ' '
    preproc += ' -DSUPERLU '

if int(lapack):
    message( green, "LAPACK support is enabled.")
    if ( arpack == 0):
        #dont repeat the message if its already been given for arpack
        message( red, "I'm assuming a GNU compiler & hard-wiring -DgFortran for")
        message( red, "cfortran to do its magic with.")
    preproc += ' -DLAPACK -DgFortran'
    libs_str += blas_lib + ' ' + lapack_lib + ' '

if int(static):
    message( yellow, "Building the static library." )
else:
    message( yellow, "Building the shared library." )
    message( yellow, "Its up to you to make sure the library is in the LD_LIBRARY_PATH.")

if int(debug):
    message( red, "DEBUG messages enabled.")
    # overwrite opts to remove any optimisation
    opts = ' '
    debug_symbols = 1
    preproc += ' -DDEBUG '

if int(paranoid):
    message( red, "PARANOID checking enabled.")
    preproc += ' -DPARANOID '

if (int(warn) > 0) :
    message( green, "WARNING level is set to 1.")
    opts += ' -Wall '

if (int(warn) >= 2) :
    message( green, "WARNING level increased to 2.")
    opts += ' -Winline '

if int(time):
    message( blue, "TIMING of selected example problems will be performed.")
    preproc += ' -DTIME '

if int(profile):
    message( blue, "PROFILING is enabled.")
    opts += ' -pg '
    link_flags += ' -pg '

if int(debug_symbols):
    message( blue, "DEBUG symbols enabled.")
    preproc += ' -g '
print

message( blue, " -----  Checking for dependencies.")

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

if int(arpack):
    arpack = fortran = 1
    if not conf.CheckLib('gfortran'):
        message( yellow, "CppNoddy build script assumes you have a fairly modern")
        message( yellow, "installation, with up-to-date GCC compiler with gfortran")
        message( red, "No libgfortran!")
        fortran = 0
    if not conf.CheckLib('arpack'):
        message( red, "No libarpack!")
        arpack = 0
    if ( arpack * fortran == 0 ):
        message( red, "ARPACK support has failed.")
        Exit(1)
    else:
        message( green, "Found ARPACK and including support for some routines.")

if int(superlu):
    # superlu => blas
    superlu = blas = 1
    if not conf.CheckLib( superlu_lib ):
        message( red, "No libsuperlu!")
        superlu = 0
    if not conf.CheckLib('blas'):
        message( red, "No libblas!")
        blas = 0
    if ( superlu * blas == 0 ):
        message( red, "SUPERLU support has failed.")
        Exit(1)
    else:
        message( green, "Found SUPERLU and including support for sparse matrix solvers.")

env = conf.Finish()

libs = libs_str.split()

message( blue, " -----  Building.")


# Build the library in ./lib
if int(static):
  env.StaticLibrary('lib/CppNoddy', src )
else:
  env.SharedLibrary('lib/CppNoddy', src, LINKFLAGS = link_flags )

# Call the Examples' SConscript if examples != 0 -- note OPTS & LFLAGS are exported to the build
if int(examples):
    SConscript('Examples/SConscript', exports='env opts preproc link_flags incdir libdir topdir libs private' )
