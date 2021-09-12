# Alchemical Transfer Method Meta-Force OpenMM Plugin
An OpenMM plugin that implements the Alchemical Transfer Potential (ATM)

[Alchemical Transfer Approach to Absolute Binding Free Energy Estimation.  J. Chem. Theory and Comput. 17, 3306 (2021)](https://doi.org/10.1021/acs.jctc.1c00266)

[Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method. ArXiv Preprint 2107.05153 (2021)](https://arxiv.org/abs/2107.05153)

## Credits

This software is written and maintained by Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>.

The plugin interface is based on the [openmmexampleplugin](https://github.com/peastman/openmmexampleplugin) by Peter Eastman.

The design of the Meta-Force is based on the implementation of OpenMM's CustomCVForce by Peter Eastman. Peter Eastman also guided much of the development of this plugin. See OpenMM's issue [#3045](https://github.com/openmm/openmm/issues/3045) for an account.

This implementation is essentially a Force-based port of the Integrator-based approach of the [openmm_sdm_plugin](https://github.com/rajatkrpal/openmm_sdm_plugin) by Rakat K. Pal and others.

Support from the National Science Foundation [CAREER 1750511](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1750511&HistoricalAwards=false) is gratefully acknowledged. 


## Requirements

OpenMM 7.5.0 or later. Tested with OpenMM 7.5.0. 

Only the Reference and OpenCL platforms are currently supported. CUDA support will be added in the near future through OpenMM's Common Compute platform.

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
mkdir build_openmm_atmplugin
cd build_openmm_atmplugin
ccmake -i ../openmm-atmmetaforce-plugin
```

Hit `c` (configure) until all variables are correctly set, then `g` to generate the makefiles. `OPENMM_DIR` should point to an existing OpenMM installation. `CMAKE_INSTALL_PREFIX` normally is the same as `OPENMM_DIR`. The plugin requires the python API. You need `python` and `swig` to install it.

Once the configuration is done do:

```
make
make install
make PythonInstall
```

The last two steps may need superuser access depending on the installation target. It is recommended to build the plugin under a `conda` environment to install the python modules without superuser access.

## Test

Go to `examples/` and edit `runopenmm` to reflect your installation and execution environment. Then:

```
./runopenmm test.py
./runopenmm test_explicit.py
```

Both tests produce a short MD run for a host-guest system at the intermediate state at λ=1/2. The first is in vacuum and the second with explicit solvent.

