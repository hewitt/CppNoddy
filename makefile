all:	
	env scons -j2  col=1 debug_symbols=1 debug=0 examples=1 time=1 superlu=1 lapack=1 paranoid=0 profile=0 private=0 static=1 warn=1

examples:
	env scons  col=0 time=1 

examples-lapack:
	env scons  col=0 time=1 lapack=1

clean:	
	env scons -c col=0

debug:
	env scons  col=0 debug=1 paranoid=1 examples=1 lapack=1

profile:
	env scons  col=0 time=1 private=1 static=1 profile=1 lapack=1

dox:
	env scons  col=0 doc=1
