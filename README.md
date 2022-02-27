# levmar

Standalone implementation in C language of the Levenberg-Marquard minimization algorithm.
No dependencies upon GSL or Eigen, nor Lapack (hence standalone).
Source code adapted and refactored from https://github.com/jturney/levmar

List of modifications:

- refactor cmake build system
- add pkg-config file
- refactor install step
- add `configure` script which allows to levmar embedded inside an autotool-based buildsytem by calling m4 macro `AX_SUBDIRS_CONFIGURE``; this configure script calls under the hood the cmake command, and fools the autotools by faking an autotool build.
- add demo applications
- add cmake option to enable internal lapack implementation (so that levmar is really self-contained, standalone).
- allow cmake to cross-compile for an embedded target (e.g. [gr712](https://www.gaisler.com/index.php/products/components/gr712rc)) :
  ```shell
  mkdir build/gr712-bare-release
  cd build/gr712-bare-release
  cmake -DCMAKE_TOOLCHAIN_FILE=../../cmake/toolchain_bcc.cmake -DCMAKE_CXX_FLAGS="-mcpu=v8 -mhard-float" -DCMAKE_BUILD_TYPE=Release -DBUILD_FPIC=OFF -DENABLE_LAPACK=OFF ../..
  make
  ```
In term of performance, it is of course less interesting than GSL or Eigen's implementation, but still interesting for embedded system use, if Eigen or GSL is not an option.
