Tortoize
========

Application to calculate ramachandran z-scores.

Building
--------

Make sure you install [libcif++](https://github.com/PDB-REDO/libcifpp) and [libzeep](https://github.com/mhekkel/libzeep) first before building.

After that, building should be as easy as typing:

```
git clone https://github.com/PBD-REDO/tortoize.git
cd tortoize
mkdir build
cd build
cmake ..
cmake --build . --config Release
ctest -C Release
cmake --install .
```

This will install the `tortoize` program in `$HOME/.local/bin`. If you want to
install elsewhere, specify the prefix with the [CMAKE_INSTALL_PREFIX](https://cmake.org/cmake/help/v3.21/variable/CMAKE_INSTALL_PREFIX.html) variable.

Usage
-----

See [manual page](doc/tortoize.pdf) for more info.
