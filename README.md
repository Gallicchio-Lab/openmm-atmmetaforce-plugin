# Alchemical Transfer Method Meta-Force OpenMM Plugin
An OpenMM plugin that implements the Alchemical Transfer Potential (ATM)

[Alchemical Transfer Approach to Absolute Binding Free Energy Estimation. arXiv:2101.07894 (2021)](https://arxiv.org/abs/2101.07894)

## Credits

This software is written and maintained by Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>.

The plugin interface is based on the [openmmexampleplugin](https://github.com/peastman/openmmexampleplugin) by Peter Eastman.

The design of the Meta-Force follows closely the implementation of OpenMM's CustomCVForce by Peter Eastman. Peter Eastman guided much of the developemnt of this plugin. See OpenMM's issue [#3045](https://github.com/openmm/openmm/issues/3045) for an account.  

Support from the National Science Foundation [CAREER 1750511](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1750511&HistoricalAwards=false) is gratefully acknowledged. 


## Requirements

OpenMM 7.5.0 or later. Tested with OpenMM 7.5.0.

## License

This software is released under the LGPL license. See LICENSE.

## Installation

Locate the OpenMM installation directory, otherwise it will default to `/usr/local/openmm`. See above regarding patching OpenMM's desmond file reader.

Download this package from github:

```
git clone https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin.git
```


Build and install the plugin with cmake. Assuming a unix system:

```
mkdir build_openmm_sdm_plugin
cd build_openmm_sdm_plugin
ccmake -i ../openmm_sdm_plugin
```

Hit `c` (configure) until all variables are correctly set, then `g` to generate the makefiles. `OPENMM_DIR` should point to an existing OpenMM installation. `CMAKE_INSTALL_PREFIX` normally is the same as `OPENMM_DIR`. The SDM plugin requires the python API. You need `python` and `swig` to install it.

Once the configuration is done do:

```
make
make install
make PythonInstall
```

The last two steps may need superuser access depending on the installation target. It is recommended to to build the plugin under a `conda` environment to install the python modules without superuser access.

## Test

Edit `runopenmm` to reflect your installation and execution environment. Then:

```
./runopenmm test.py
./runopenmm test_explicit.py
```

Both tests produce a short MD run for a host-guest system at the intermediate state at λ=1/2. The first is in vacuum and the second with explicit solvent.

