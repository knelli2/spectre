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
  module load gmp/6.2.1-gcc-13.2.0-lcnhyse
  module load mpfr/4.2.0-gcc-13.2.0-yy2fkq5
  module load mpc/1.3.1-gcc-13.2.0-5kgoftq
  module load zlib-ng/2.1.3-gcc-13.2.0-jetnfwa
  module load zstd/1.5.5-gcc-13.2.0-t2lua3l
#  module load llvm/14.0.6-gcc-11.3.1-3o7col4
  module load gcc/13.2.0-gcc-13.2.0-w55nxkl
  module load cmake/3.20.2-gcc-13.2.0-rp74vpv
  module load libfabric/1.18.1-gcc-11.3.1-3xjzfrf
  module load numactl/2.0.14-gcc-11.3.1-dztkbhb
  module load openssh/8.7p1-gcc-11.3.1-b4y76zs
  module load slurm/22.05.6-gcc-11.3.1-atym2np
  module load zlib-ng/2.1.3-gcc-11.3.1-77c4j3n
  module load openmpi/4.1.5-gcc-11.3.1-tzth463
  module load pkgconf/1.7.3-gcc-11.3.1-x2l7dmy
  module load hdf5/1.14.2-gcc-11.3.1-olj2klq
  module load intel-oneapi-mpi/2021.10.0-gcc-11.3.1-xj5ixri
  module load libpciaccess/0.17-gcc-13.2.0-r2cijnn
  module load libiconv/1.17-gcc-13.2.0-ntov4te
  module load xz/5.4.1-gcc-13.2.0-4xkm5xo
  module load libxml2/2.10.3-gcc-13.2.0-fr6jcjz
  module load ncurses/6.4-gcc-13.2.0-4o2yj6n
  module load hwloc/2.9.1-gcc-13.2.0-gzvfolk
  module load intel-tbb/2021.9.0-gcc-13.2.0-6nwk3ml
  module load intel-oneapi-mkl/2023.2.0-gcc-13.2.0-ohvyk7g
  module load bzip2/1.0.8-gcc-13.2.0-ypy2y5b
  module load libmd/1.0.4-gcc-13.2.0-s3xkdtu
  module load libbsd/0.11.7-gcc-13.2.0-cl7lha3
  module load expat/2.5.0-gcc-13.2.0-hymnqig
  module load readline/8.2-gcc-13.2.0-pinqxvi
  module load gdbm/1.23-gcc-13.2.0-exyzfc5
  module load tar/1.34-gcc-13.2.0-3nxpzly
  module load gettext/0.21.1-gcc-13.2.0-yiohc4s
  module load libffi/3.4.4-gcc-13.2.0-6hw6v3w
  module load libxcrypt/4.4.35-gcc-13.2.0-4fcttz6
  module load openssl/3.1.2-gcc-13.2.0-egeilil
  module load sqlite/3.42.0-gcc-13.2.0-wri5a47
  module load util-linux-uuid/2.38.1-gcc-13.2.0-2dkujln
  module load python/3.10.12-gcc-13.2.0-xhirdjb
  module load yaml-cpp/0.7.0-gcc-13.2.0-szs4phw
  module load gsl/2.7.1-gcc-13.2.0-6hld2nf
  module use /central/groups/sxs/knelli/tools/amd/modules/
  module load blaze/3.8
  module load boost/1.82.0
  module load brigand/master
  module load catch/3.4.0
  module load charm/7.0.0
  module load libsharp/1.0.0
  module load libxsmm/1.16.1
  module load spectre-python/2023-12
}

spectre_unload_modules() {
  module unload gmp/6.2.1-gcc-13.2.0-lcnhyse
  module unload mpfr/4.2.0-gcc-13.2.0-yy2fkq5
  module unload mpc/1.3.1-gcc-13.2.0-5kgoftq
  module unload zlib-ng/2.1.3-gcc-13.2.0-jetnfwa
  module unload zstd/1.5.5-gcc-13.2.0-t2lua3l
 # module unload llvm/14.0.6-gcc-11.3.1-3o7col4
  module unload gcc/13.2.0-gcc-13.2.0-w55nxkl
  module unload cmake/3.20.2-gcc-13.2.0-rp74vpv
  module unload libfabric/1.18.1-gcc-11.3.1-3xjzfrf
  module unload numactl/2.0.14-gcc-11.3.1-dztkbhb
  module unload openssh/8.7p1-gcc-11.3.1-b4y76zs
  module unload slurm/22.05.6-gcc-11.3.1-atym2np
  module unload zlib-ng/2.1.3-gcc-11.3.1-77c4j3n
  module unload openmpi/4.1.5-gcc-11.3.1-tzth463
  module unload pkgconf/1.7.3-gcc-11.3.1-x2l7dmy
  module unload hdf5/1.14.2-gcc-11.3.1-olj2klq
  module unload intel-oneapi-mpi/2021.10.0-gcc-11.3.1-xj5ixri
  module unload libpciaccess/0.17-gcc-13.2.0-r2cijnn
  module unload libiconv/1.17-gcc-13.2.0-ntov4te
  module unload xz/5.4.1-gcc-13.2.0-4xkm5xo
  module unload libxml2/2.10.3-gcc-13.2.0-fr6jcjz
  module unload ncurses/6.4-gcc-13.2.0-4o2yj6n
  module unload hwloc/2.9.1-gcc-13.2.0-gzvfolk
  module unload intel-tbb/2021.9.0-gcc-13.2.0-6nwk3ml
  module unload intel-oneapi-mkl/2023.2.0-gcc-13.2.0-ohvyk7g
  module unload bzip2/1.0.8-gcc-13.2.0-ypy2y5b
  module unload libmd/1.0.4-gcc-13.2.0-s3xkdtu
  module unload libbsd/0.11.7-gcc-13.2.0-cl7lha3
  module unload expat/2.5.0-gcc-13.2.0-hymnqig
  module unload readline/8.2-gcc-13.2.0-pinqxvi
  module unload gdbm/1.23-gcc-13.2.0-exyzfc5
  module unload tar/1.34-gcc-13.2.0-3nxpzly
  module unload gettext/0.21.1-gcc-13.2.0-yiohc4s
  module unload libffi/3.4.4-gcc-13.2.0-6hw6v3w
  module unload libxcrypt/4.4.35-gcc-13.2.0-4fcttz6
  module unload openssl/3.1.2-gcc-13.2.0-egeilil
  module unload sqlite/3.42.0-gcc-13.2.0-wri5a47
  module unload util-linux-uuid/2.38.1-gcc-13.2.0-2dkujln
  module unload python/3.10.12-gcc-13.2.0-xhirdjb
  module unload yaml-cpp/0.7.0-gcc-13.2.0-szs4phw
  module unload gsl/2.7.1-gcc-13.2.0-6hld2nf
  module unload blaze/3.8
  module unload boost/1.82.0
  module unload brigand/master
  module unload catch/3.4.0
  module unload charm/7.0.0
  module unload libsharp/1.0.0
  module unload libxsmm/1.16.1
  module unload spectre-python/2023-12
  module unuse /central/groups/sxs/knelli/tools/amd/modules/
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
          -D MEMORY_ALLOCATOR=SYSTEM \
          -D BUILD_PYTHON_BINDINGS=ON \
          -D MACHINE=CaltechHpc \
          "$@" \
          $SPECTRE_HOME
}
