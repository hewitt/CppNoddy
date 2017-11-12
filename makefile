all:
	env scons -j3  col=1 debug=0 examples=1 time=0 petsc=1 lapack=1 slepc=0 paranoid=0 profile=0 private=0

clean:
	env scons -c

dox:
	env scons  col=0 doc=1
