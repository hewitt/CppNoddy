touch ./validate.log
rm -rf ./validate.log
 cd .././Arclength
.././Arclength/Arc_Shoot_FalknerSkan.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Arclength
.././Arclength/Arc_Trans_Fold.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Arclength
.././Arclength/Arc_circle.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Arclength
.././Arclength/Arc_circle_vector.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Berman.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Blasius.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Harmonic.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Karman.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Karman_Jacobian.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Karman_adapt.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Karman_arc.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_Troesch.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/BVP_nonIdentity.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/Shoot_Berman.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././BVP
.././BVP/Shoot_Blasius.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/MPI_MUMPS.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/MatrixMult.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/MatrixSolves.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/MatrixSparseSolves.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/Matrix_CompareNative.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/Vec_DenseVector.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/Vec_Overloading.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Containers
.././Containers/Vec_SparseVector.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_Harmonic.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_Harmonic_easy.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_Harmonic_sparse.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_OrrSommerfeld.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_OrrSommerfeld_easy.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_OrrSommerfeld_neutralcurve.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_OrrSommerfeld_sparse.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_OrrSommerfeld_student.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_Rayleigh.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_Shoot_Biharmonic.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_complex.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_complex_sparse.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_local_Harmonic.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././EVP
.././EVP/EVP_real_sparse.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Generic
.././Generic/ExceptionChecks.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Generic
.././Generic/NewtonIter.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Generic
.././Generic/Quad.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Generic
.././Generic/TrivialComplex.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_1D
.././HYP_1D/HYP_acoustic_impedance.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_1D
.././HYP_1D/HYP_acoustic_reflection.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_1D
.././HYP_1D/HYP_nonlinear_advection.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_1D
.././HYP_1D/HYP_radial_dam_break.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_1D
.././HYP_1D/HYP_shocktube_Sod.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_linear_advection_xy.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_nonlinear_advection_x.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_nonlinear_advection_xy.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_nonlinear_advection_y.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_radial_dam_break.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././HYP_2D
.././HYP_2D/HYP_2D_shallow_source.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IBVP
.././IBVP/IBVP_Karman.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IBVP
.././IBVP/IBVP_diffusion.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IBVP
.././IBVP/IBVP_diffusion_nonlinear.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IBVP
.././IBVP/IBVP_nonlinear_advdiff.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IVP
.././IVP/IVP_Harmonic.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././IVP
.././IVP/IVP_Lorenz.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Meshes
.././Meshes/1DNodeMesh.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Meshes
.././Meshes/1DNodeMesh_Airy.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Meshes
.././Meshes/2DNodeMesh.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Poisson
.././Poisson/Poisson_C.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Poisson
.././Poisson/Poisson_Stokes.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././Poisson
.././Poisson/Poisson_m.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././double_IBVP
.././double_IBVP/IBVP_linear.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././double_IBVP
.././double_IBVP/IBVP_nonlinear_fast.out 2>> error.log | tee -a ../Validation/validate.log
 cd .././double_IBVP
.././double_IBVP/IBVP_nonlinear_slow.out 2>> error.log | tee -a ../Validation/validate.log
cd ../Validation 
echo 
echo '=================== TEST RESULTS ===================='
echo 
a=`grep tee ./validate.sh | wc -l`
echo -n 'Number of tests run              = '
echo $(($a-1))
echo -n 'Number of tests returning a pass = '
grep PASSED ./validate.log | wc -l 
echo -n 'Number of tests returning a skip = '
grep SKIPPED ./validate.log | wc -l 
echo -n 'Number of tests returning a fail = '
grep FAILED ./validate.log | wc -l 
