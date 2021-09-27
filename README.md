# Alchemical Transfer Method Meta-Force OpenMM Plugin
An OpenMM plugin that implements the Alchemical Transfer Potential (ATM) for absolute and relative binding free energy calculations. 

## References

The Alchemical Transfer Method (ATM) is based on a coordinate transformation to transfer directly the ligand from the solvent bulk to the binding site. A similar coordinate transformation is used to swap the positions of two ligands between the solvent and the binding site for relative binding free energy calculations. The method uses one simulation box, prepared with standard tools, without the need of dummy atoms nor alchemical topologies. In principle, any molecular file format supported by OpenMM is supported by this plugin. ATM and its applications are described in the following publications:

[Alchemical Transfer Approach to Absolute Binding Free Energy Estimation.  J. Chem. Theory and Comput. 17, 3306 (2021)](https://doi.org/10.1021/acs.jctc.1c00266)

[Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method. ArXiv Preprint 2107.05153 (2021)](https://arxiv.org/abs/2107.05153)

[Application of the Alchemical Transfer and Potential of Mean Force Methods to the SAMPL8 Host-Guest Blinded Challenge](https://arxiv.org/abs/2107.05155)

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

Go to `example/abfe/`, set the `OPENMM_DIR` environment variable to point to the OpenMM installation location with this plugin. Then do:

```
../runopenmm abfe.py
```

Repeat for the example under `example/rbfe/`:

```
../runopenmm rbfe.py
```

Each test produces a short MD run for a [SAMPL8 host-guest system](https://arxiv.org/abs/2107.05155) at the alchemical intermediate state at λ=1/2. The first is for the absolute binding free energy calculation and the second for relative binding. Each generates an output file, `.out`, of the form:
```
00.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69967.037957 -0.956023
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69923.338273 -0.836520
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69737.383816 5.048996
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69704.356182 1.195029
...
```
where each line is a sample from the MD run. The first column is the temperature, the second is λ, then λ1,  λ2, α, u0, w0, the total potential energy, and the last column is the alchemical perturbation energy u in kcal/mol energy units. The set of perturbation energies collected at a series of λ-values is used to compute the binding free energy. See the references above for the theoretical background and the details of the method.  
