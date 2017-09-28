cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(votca_xtp_gpu LANGUAGES CXX CUDA)

file(GLOB_RECURSE VOTCA_GPU_SOURCES *.cu)

add_library(votca_xtp_gpu OBJECT
  ${VOTCA_GPU_SOURCES}
  ../tools/GPUNumericalIntegrationTest.cu
  )

target_compile_features(votca_xtp_gpu PUBLIC cxx_std_11)
set_property(TARGET votca_xtp_gpu PROPERTY CUDA_STANDARD 11)

set_target_properties(votca_xtp_gpu PROPERTIES  POSITION_INDEPENDENT_CODE ON)
