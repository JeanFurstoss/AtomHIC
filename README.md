### AtomHic library : generate your own executable to post-treat or generate atomic systems
- Required
   - CMake 2.8.9 or higher
   - C++ compiler: gcc 8 or higher 

- Compilation
  - Linux
      - mkdir build
      - cd build
      - cmake ..
      - ccmake .. (generate shared lib and tools) 
      - make
- Run
   - Linux
     - Launch "build/bin/exe.." executable

### No available yet

- Tests
   - Linux
     - cmake .. -DBUILD_TESTS=true
     - ctest -V

### some documentation for using the library and generating executable is present in doc/

Some part of the code are parallelized using openmpi, to benefit for this, use : export OMP_NUM_THREADS="number of thread you want to use"

Warning, for the moment the default boundary conditions are periodic, for systems with non periodic boundaries the better is to increase the cell size to dimensions higher than cutoff radius TODO
