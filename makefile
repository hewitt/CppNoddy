all:
	env scons -j1  col=1 debug_symbols=1 debug=0 examples=1 time=0 superlu=0 mumps=1 petsc=1 lapack=1 slepc=1 paranoid=0 profile=0 private=1

clean:
	env scons -c

dox:
	env scons  col=0 doc=1
