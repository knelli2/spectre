#!/bin/env sh

# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_setup_modules() {
    echo "All modules on Wheeler are provided by the system"
}

spectre_unload_modules() {
    module unload gcc/7.3.0
    module unload blaze/3.8
    module unload boost/1.65.0-gcc-6.4.0
    module unload brigand/master
    module unload catch/2.13.3
    module unload gsl/2.1
    module unload libsharp/1.0.0
    module unload libxsmm/1.16.1
    module unload openblas/0.2.18
    module unload papi/5.5.1
    module unload yaml-cpp/master
    module unload impi/2017.1
    module unload cmake/3.18.2
    module unload ninja/1.10.0
    module unload doxygen/1.8.13
    module unload git/2.8.4
    module unload llvm/10.0.0
    module unload charm/6.10.2-intelmpi-smp
    module unload python/anaconda3-2019.10
    module unload pybind11/2.6.1
}

spectre_load_modules() {
    module load gcc/7.3.0
    module load blaze/3.8
    module load boost/1.65.0-gcc-6.4.0
    module load brigand/master
    module load catch/2.13.3
    module load gsl/2.1
    module load libsharp/1.0.0
    module load libxsmm/1.16.1
    module load openblas/0.2.18
    module load papi/5.5.1
    module load yaml-cpp/master
    module load impi/2017.1
    module load cmake/3.18.2
    module load ninja/1.10.0
    module load doxygen/1.8.13
    module load git/2.8.4
    module load llvm/10.0.0
    module load charm/6.10.2-intelmpi-smp
    module load python/anaconda3-2019.10
    module load pybind11/2.6.1
}

spectre_run_cmake() {
    if [ -z ${SPECTRE_HOME} ]; then
        echo "You must set SPECTRE_HOME to the cloned SpECTRE directory"
        return 1
    fi
    spectre_load_modules
    # Notes:
    # - Set CMAKE_PREFIX_PATH to pick up packages consistent with the anaconda
    #   module, such as zlib. The anaconda module on Wheeler does not set this
    #   automatically.
    cmake -D CHARM_ROOT=$CHARM_ROOT \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_C_COMPILER=clang \
          -D CMAKE_CXX_COMPILER=clang++ \
          -D CMAKE_Fortran_COMPILER=gfortran \
          -D MEMORY_ALLOCATOR=SYSTEM \
          -D BUILD_PYTHON_BINDINGS=OFF \
          -D CMAKE_PREFIX_PATH="$PYTHON_HOME" \
          "$@" \
          $SPECTRE_HOME
}

spectre_run_cmake_charm7() {
    if [ -z ${SPECTRE_HOME} ]; then
        echo "You must set SPECTRE_HOME to the cloned SpECTRE directory"
        return 1
    fi
    spectre_load_modules
    CHARM_ROOT="/home/knelli/tools/charm_7/mpi-linux-x86_64-smp"
    # Notes:
    # - Set CMAKE_PREFIX_PATH to pick up packages consistent with the anaconda
    #   module, such as zlib. The anaconda module on Wheeler does not set this
    #   automatically.
    cmake -D CHARM_ROOT=$CHARM_ROOT \
          -D CMAKE_BUILD_TYPE=Release \
          -D CMAKE_C_COMPILER=clang \
          -D CMAKE_CXX_COMPILER=clang++ \
          -D CMAKE_Fortran_COMPILER=gfortran \
          -D MEMORY_ALLOCATOR=SYSTEM \
          -D BUILD_PYTHON_BINDINGS=on \
          -D CMAKE_PREFIX_PATH="$PYTHON_HOME" \
          -D ASAN=ON \
          "$@" \
          $SPECTRE_HOME
}
