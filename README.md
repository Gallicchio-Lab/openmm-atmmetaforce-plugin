#### {#mainpage}
# Alchemical Transfer Method Meta-Force OpenMM Plugin

An OpenMM plugin that implements the Alchemical Transfer Potential (ATM) for absolute and relative binding free energy calculations. 

## References

The Alchemical Transfer Method (ATM) is based on a coordinate transformation to transfer directly the ligand from the solvent bulk to the binding site. A similar coordinate transformation is used to swap the positions of two ligands between the solvent and the binding site for relative binding free energy calculations. The method uses one simulation box, prepared with standard tools, without the need of dummy atoms nor alchemical topologies. In principle, any molecular file format supported by OpenMM is supported by this plugin. 

ATM and its implementation are described in the open access article:

Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.  [J. Chem. Inf. Model.  62, 309 (2022)](https://doi.org/10.1021/acs.jcim.1c01129)

Refer to the publication above for a detailed description of the ATM method and the parameters used in this API 
and please cite it to support our work if you use this software in your research.

## Credits

This software is written and maintained by Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>.

The plugin interface is based on the [openmmexampleplugin](https://github.com/peastman/openmmexampleplugin) by Peter Eastman.

The design of the Meta-Force is based on the implementation of OpenMM's CustomCVForce by Peter Eastman. Peter Eastman also guided much of the development of this plugin. See OpenMM's issue [#3045](https://github.com/openmm/openmm/issues/3045) for an account.

This implementation is essentially a Force-based port of the Integrator-based approach of the [openmm_sdm_plugin](https://github.com/rajatkrpal/openmm_sdm_plugin) by Rakat K. Pal and others.

Support from the National Science Foundation [CAREER 1750511](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1750511&HistoricalAwards=false) is gratefully acknowledged. 


## Requirements

OpenMM 7.7. Releases prior to 0.3.1 support 7.5.0 and 7.6.0.

## License

This software is released under the LGPL license. See LICENSE.

## Installation

### From conda-forge

The latest stable binary release is available in conda-forge. To install it do:

```
conda install -c conda-forge openmm-atmmetaforce-plugin
```

### From sources

Locate the OpenMM installation directory, otherwise it will default to `/usr/local/openmm`.

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

## Documentation

The ATM Meta Force plugin has a C++ API that is documented inline in the `openmmapi/include/ATMMetaForce.h` header file.
The Python API mirrors the C++ API and is accessed by importing the `atmmetaforce.py` module (see example below).
The Python interface also includes a `ATMMetaForceUtils` module
that includes utility functions, such as restraining potentials, often used in alchemical binding calculations.

The latest documentation of the C++ API and the `ATMMetaForceUtils` Python module can be viewed in html and latex using `doxygen`
```
conda install -c conda-forge doxygen
```
Download the source from this repository, then, for example,
```
cd openmm-atmmetaforce-plugin/doc
doxygen
google-chrome html/index.html
```

A [pdf](https://drive.google.com/file/d/1FFDYBEjTh6VUWg9YnDZzMOOhZLj_dr6L/view?usp=sharing) version of the API documentation is available.

## Example of Usage (v0.3.0)

Rather than implementing a standard molecular interaction potential, the ATMMetaForce class provides
an infrastructure to compute a chosen set of OpenMM Forces before and after the application of the
transformation of the system coordinates, and to merge the resulting energies and atomic forces
according to an alchemical potential. In this respect it is similar in spirit to OpenMM's CustomCVForce.

For example, when using a linear alchemical potential, if `U1(x)` and `U2(x)` are the potential energies from the
Forces managed by ATMMetaForce before and after the application of the coordinate transformation, the ATMMetaForce
adds to the OpenMM's system Forces the potential
```
W(x; lambda) = (1 - lambda) U1(x) + lambda U2(x)
```
where `0 < lambda < 1` is the alchemical parameter. This is for a positive alchemical direction (`direction` = 1).
If the direction of the transformation is reversed (`direction` = -1)
the alchemical potential becomes `U(x) = (1 - lambda) U2(x) + lambda U1(x)` with `U1(x)` and `U2(x)` reversed.

In practice, 
a more elaborate softplus hybrid alchemical potential is used (see the paper above) which depends on two lambda parameters,
`lambda1` and `lambda2`, as well as on other parameters. The softplus alchemical potential has the expression
```
W(x; lambda1, lambda2, alpha, u0, w0) =
(lambda2 - lambda1) ln{ 1 + exp[-alpha (u - u)] }/alpha + lambda2 u + w0
```
where `u = u(x) = U2(x) - U1(x)` is the perturbation energy (or `U1(x) - U2(x)` if `direction` = -1).
Note that the softplus alchemical potential is the same
as the linear potential for `lambda2 = lambda1 = lambda`. Also note that the perturbation energy is further
modified by a softcore mapping function that caps its maximum value according to the parameters `umax`, `ubcore`, and `acore` below.

To use the class, place the Forces that should be recalculated before and after the transformation
in one or more force groups. Create an ATMMetaForce Force object, assign an unused Force group
id to it, then add all of the particles to it with their displacement transformation. Then perform 
minimization/MD as usual telling the Integrator to consider only the ATMMetaForce and the Forces not 
managed by the ATMMetaForce.

For example, the following excerpt (see the `examples`) sets up a transformation in which a a ligand
(atoms 0, 1, 2, 3, 4) is displaced by 22 Angstroms in the each of the x, y, z directions and simulates
at the midpoint of the transformation (lambda = 1/2). The non-bonded force
is placed in Force group 1 and is recalculated before and after the transformation. Forces, such as as the bonded Forces,
placed in default Forge group 0 are not affected. The ATM Meta Force is assigned to group 2. The integrator is then told
to consider only the bonded Forces and the ATM Meta Force (groups 0 and 2) when integrating the equations of motion/minimization etc. The non-bonded force is skipped because it is managed by the ATM Meta Force.

```
from atmmetaforce import *

atm_utils = ATMMetaForceUtils(system)

lmbd = 0.5
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole
umsc =  200.0 * kilocalorie_per_mole
ubcore = 100.0 * kilocalorie_per_mole
acore = 0.062500
direction = 1.0

atmforcegroup = 2
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)

atmvariableforcegroups = [nonbonded_force_group]
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, direction, atmvariableforcegroups )

for at in top.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
lig_atoms = [0, 1, 2, 3, 4]
displ = [22.0, 22.0, 22.0]
for i in lig_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
atmforce.setForceGroup(atmforcegroup)

system.addForce(atmforce)
print("Using ATM Meta Force plugin version = %s" % ATMMETAFORCE_VERSION)

integrator = LangevinIntegrator(temperature, frictionCoeff, MDstepsize)
integrator.setIntegrationForceGroups({0,atmforcegroup})

```

## Tests & Further Examples

Download the source distribution (see above), go to `example/abfe/`, set the `OPENMM_DIR` environment variable to point to the OpenMM installation location with this plugin. Then do:

```
../runopenmm abfe.py
```

Repeat for the example under `example/rbfe/`:

```
../runopenmm rbfe.py
```

Each test produces a short MD run for a [SAMPL8 host-guest system](https://arxiv.org/abs/2107.05155) at the alchemical intermediate state at lambda=1/2. The first is for the absolute binding free energy calculation and the second for relative binding. Each generates an output file, `.out`, of the form:
```
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69967.037957 -0.956023
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69923.338273 -0.836520
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69737.383816 5.048996
300.000000 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 -69704.356182 1.195029
...
```
where each line is a sample from the MD run. The first column is the temperature, the second (see below) is lambda, then lambda1,  lambda2, alpha, u0, w0, the total potential energy, and the last column is the alchemical perturbation energy u in kcal/mol energy units. The set of perturbation energies collected at a series of lambda-values is used to compute the binding free energy. See the references above for the theoretical background and the details of the method.  

