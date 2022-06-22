c = run_command('ls_all.sh')
sources = c.stdout().strip().split('\n')

# Build with PETSC hooks
if slepc
   CppNoddy = shared_library('CppNoddy',
	 sources,
	 include_directories : [ inc ],
	 dependencies : [ lapack_dep, petsc_dep, slepc_dep, mpich_dep ],
	 cpp_args : build_args,
	 install : true ) 
elif petscd or petscz
   CppNoddy = shared_library('CppNoddy',
	 sources,
	 include_directories : [ inc ],
	 dependencies : [ lapack_dep, petsc_dep, mpich_dep ],
	 cpp_args : build_args,
	 install : true )  
else
   CppNoddy = shared_library('CppNoddy',
	 sources,
	 include_directories : [ inc ],
	 dependencies : [ lapack_dep ],
	 cpp_args : build_args,
	 install : true )
endif

