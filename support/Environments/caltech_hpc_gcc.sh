#!/bin/env sh

# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_setup_modules() {
    module use /central/groups/sxs/modules
    echo "Place the following line in you '~/.bashrc' so you don't have to "
    echo "run 'spectre_setup_modules' every time you log in:"
    echo ""
    echo "module use /central/groups/sxs/modules"
}

spectre_load_modules() {
    # System-installed
    module use /central/groups/sxs/knelli/modules/
    module load git/2.37.2
    module load compilers/binutils/2.41
    module load compilers/gcc/11.4.0
    module load tools/cmake/3.28.1
    module load libraries/blaze/3.8
    module load libraries/boost/1.82.0
    module load libraries/brigand/master
    module load libraries/catch/3.4.0
    module load oneapi/mpi-2021.11
    module load libraries/charm/7.0.0
    module load libraries/gsl/2.7
    module load libraries/hdf5/1.12.3
    module load libraries/jemalloc/5.3.0
    module load libraries/libsharp/1.0.0
    module load libraries/libxsmm/1.16.1-avx2
    module load libraries/openblas/0.3.25-avx2
    # module load libraries/libxsmm/1.16.1-avx512
    # module load libraries/openblas/0.3.25-avx512
    module load libraries/yaml-cpp/0.8.0
    module load languages/python/spectre-python
}

spectre_unload_modules() {
    module use /central/groups/sxs/knelli/modules/
    module unload git/2.37.2
    module unload compilers/binutils/2.41
    module unload compilers/gcc/11.4.0
    module unload tools/cmake/3.28.1
    module unload libraries/blaze/3.8
    module unload libraries/boost/1.82.0
    module unload libraries/brigand/master
    module unload libraries/catch/3.4.0
    module unload oneapi/mpi-2021.11
    module unload libraries/charm/7.0.0
    module unload libraries/gsl/2.7
    module unload libraries/hdf5/1.12.3
    module unload libraries/jemalloc/5.3.0
    module unload libraries/libsharp/1.0.0
    module unload libraries/libxsmm/1.16.1-avx2
    module unload libraries/openblas/0.3.25-avx2
    # module unload libraries/libxsmm/1.16.1-avx512
    # module unload libraries/openblas/0.3.25-avx512
    module unload libraries/yaml-cpp/0.8.0
    module unload languages/python/spectre-python
}

spectre_run_cmake() {
    if [ -z ${SPECTRE_HOME} ]; then
        echo "You must set SPECTRE_HOME to the cloned SpECTRE directory"
        return 1
    fi
    spectre_load_modules
    cmake -D CMAKE_C_COMPILER=gcc \
          -D CMAKE_CXX_COMPILER=g++ \
          -D CMAKE_Fortran_COMPILER=gfortran \
          -D CHARM_ROOT=$CHARM_ROOT \
          -D CMAKE_BUILD_TYPE=Release \
          -D MEMORY_ALLOCATOR=JEMALLOC \
          -D BUILD_PYTHON_BINDINGS=ON \
          -D MACHINE=CaltechHpc \
          -D OVERRIDE_ARCH=skylake \
          -D YAMLCPP_ROOT=/central/groups/sxs/knelli/libraries/yaml-cpp/0.8.0 \
          "$@" \
          $SPECTRE_HOME
}
