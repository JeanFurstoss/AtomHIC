Advanced example using c++, cmake, ctest, doxygen

- Required
   - CMake 2.8.9 or higher
   - C++ compiler (gnu, visual studio...)

- Compilation
  - Linux
      - mkdir build
      - cd build
      - cmake ..
      - make
- Run
   - Linux
     - Launch "build/bin/exe.." executable
- Tests
   - Linux
     - cmake .. -DBUILD_TESTS=true
     - ctest -V
   - Windows
     - launch cmake-gui, BUILD_TESTS activated
     - generate "RUN_TESTS_VERBOSE" solution under visual

- Source documentation
   - Linux
     - cmake .. -DBUILD_DOCUMENTATION=true
     - make doc
     - firefox build/doc/html/index.html
   - Windows
     - launch cmake-gui, BUILD_DOCUMENTATION activated
     - generate "doc" solution under visual
     - firefox build/doc/html/index.html
