Tortoize
========

[![github CI](https://github.com/PDB-REDO/tortoize/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/PDB-REDO/tortoize/actions)

Application to calculate ramachandran z-scores.

Building
--------

The easiest way to install tortoize is by installing [CCP4](https://www.ccp4.ac.uk/download/index.php)

It is possible to install tortoize on Linux without having CCP4. In that case you will have install some dependencies first. On Debian this boils down to:

```console
sudo apt-get update && sudo apt-get install libcatch2-dev libeigen3-dev
```

And on Ubuntu, slightly different:

```console
sudo apt-get update && sudo apt-get install catch2 libeigen3-dev
```

After that, building and installing should be as simple as:

```console
git clone https://github.com/PDB-REDO/tortoize.git
cd tortoize
cmake -S . -B build
cmake --build build
sudo cmake --install build
```

Usage
-----

See [manual page](doc/tortoize.pdf) for more info.
