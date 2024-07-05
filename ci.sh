#!/usr/bin/env bash
# Modes
# 1: Standalone instalation outside current folder, including a fresh petsc build
# 2: Fast CI with current source code and current petsc installation
mode=2

function test(){
  echo "***************** TEST: $1 ****************** "
  local test_name=$1
  local solution_name=$2
  cd $3
  make
  ./$1.exe ${MAD_ROOT}/ci/${test_name}/par &> ${MAD_ROOT}/ci/${test_name}/log_$test_name.txt
  ground_truth=${MAD_ROOT}/ci/${test_name}/
  candidate=${MAD_RESULTS}/ci/${test_name}
  if  cmp -s $ground_truth/${solution_name} $candidate/${solution_name}; then
    return_status="SUCCES"
    printf "***************** $test_name: $return_status ***************** \n\n"
  else
    return_status="FAIL"
    printf "***************** $test_name: $return_status ***************** \n\n"
  fi

}

if [ $mode -eq 1 ]
then
  cd ci    
  git clone -b release https://gitlab.com/petsc/petsc.git ./petsc/
elif [ $mode -eq 2 ]
then
  echo "***************** BUILDING MAD ***************** "
  make clean
  make
  rm -r ${MAD_RESULTS}/ci/

  test "laplacian_square" "u.scl" "${MAD_ROOT}/projects/pdes/00_laplacian/00_square/"
  laplacian_output=$return_status
  test "heat_square" "temperature.00000.scl" "${MAD_ROOT}/projects/pdes/01_heat/00_square/"
  heat_square_output=$return_status
  test "normals_split_2D" "n_x.scl" "${MAD_ROOT}/projects/post-proc/00_normals_split/"
  nx_split_2d_output=$return_status
  test "normals_split_2D" "n_y.scl" "${MAD_ROOT}/projects/post-proc/00_normals_split/"
  ny_split_2d_output=$return_status
  test "normals_contiguous_2D" "normals.vct" "${MAD_ROOT}/projects/post-proc/01_normals_contiguous/"
  normals_contiguous_2D=$return_status
  test "normals_split_3D" "n_x.scl" "${MAD_ROOT}/projects/post-proc/02_normals_split_3D/"
  nx_split_3d_output=$return_status
  test "normals_split_3D" "n_y.scl" "${MAD_ROOT}/projects/post-proc/02_normals_split_3D/"
  ny_split_3d_output=$return_status
  test "normals_split_3D" "n_z.scl" "${MAD_ROOT}/projects/post-proc/02_normals_split_3D/"
  nz_split_3d_output=$return_status
  test "normals_contiguous_3D" "normals.vct" "${MAD_ROOT}/projects/post-proc/03_normals_contiguous_3D/"
  normals_contiguous_3D_status=$return_status
  test "von_karman" "ux.00000.scl" "${MAD_ROOT}/projects/pdes/13_navier_stokes_nonlinear/00_fpc/00_fastVonKarman/"
  von_karman_output=$return_status
  test "ns_rectangle_nn" "ux.00019.scl" "${MAD_ROOT}/projects/pdes/13_navier_stokes_nonlinear/00_fpc/02_nnFPC/"
  ns_rectangle_nn_status=$return_status
  test "ns_tube_simvascular" "wkssl_pressure.txt" "${MAD_ROOT}/projects/pdes/13_navier_stokes_nonlinear/07_tube_simvascular/"
  ns_tube_simvascular=$return_status
  test "hammer" "ctrlP.txt" "${MAD_ROOT}/projects/pdes/06_hammer/01_picard/"
  hammer_status=$return_status
  test "elasticity" "disp.vct" "${MAD_ROOT}/projects/pdes/02_elasticity/01_square/"
  elasticity_status=$return_status
  test "biot_frequency" "phi_re.scl" "${MAD_ROOT}/projects/pdes/15_complex_poro/00_square/"
  biot_frequency_status=$return_status
fi

printf "= = = = = = = = = = = = = = = = = = = =\n"
printf "\n CONTINUOUS INTEGRATION SUMMARY\n"
printf "    Test Name        |    Output\n"
printf " laplacian_square    |    $laplacian_output \n"
printf "  heat_square        |    $heat_square_output \n"
printf " normals split 2D    |    x: $nx_split_2d_output, y: $ny_split_2d_output \n"
printf " normals cont. 2D    |    $normals_contiguous_2D \n"
printf " normals split 3D    |    x: $nx_split_3d_output, y: $ny_split_3d_output, z: $nz_split_3d_output \n"
printf " normals cont. 3D    |    $normals_contiguous_3D_status \n"
printf "   von_karman        |    $von_karman_output \n"
printf " ns_rectangle_nn     |    $ns_rectangle_nn_status \n"
printf " ns_tube_simvascular |    $ns_tube_simvascular \n"
printf "       hammer        |    $hammer_status\n"
printf "     elasticity      |    $elasticity_status\n"
printf "   biot_frequency    |    $biot_frequency_status\n"
printf "= = = = = = = = = = = = = = = = = = = =\n\n"
