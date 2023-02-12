import openmm as mm
import openmm.app as mmapp
import openmm.unit as unit
import numpy as np
import pytest
import atmmetaforce as atm
from sys import stdout
import os, re, time, shutil, math
from datetime import datetime

jobname = "temoa-g1"
platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]

@pytest.mark.parametrize('prmtop_file',
                         [jobname + '.prmtop']) 
@pytest.mark.parametrize('inpcrd_file,',
                         [jobname + '.inpcrd']) 
@pytest.mark.parametrize('xml_file',
                         [jobname + '-equil.xml'])
@pytest.mark.parametrize('platform', platforms)
def test_BindingEnergy(prmtop_file, inpcrd_file, xml_file, platform):
    lmbd = 0.5
    lambda1 = lmbd
    lambda2 = lmbd
    alpha = 0.0 / unit.kilocalorie_per_mole
    u0 = 0.0 * unit.kilocalorie_per_mole
    w0coeff = 0.0 * unit.kilocalorie_per_mole
    umsc =  200.0 * unit.kilocalorie_per_mole
    ubcore = 100.0 * unit.kilocalorie_per_mole
    acore = 0.062500
    direction = 1.0

    rcpt_resid = 1
    lig_resid = 2
    displ = [ 22.0, 22.0, 22.0 ]

    displacement      = [  displ[i] for i in range(3) ] * unit.angstrom
    lig_restr_offset  = [    0      for i in range(3) ] * unit.angstrom

    prmtop = mmapp.AmberPrmtopFile(prmtop_file)
    inpcrd = mmapp.AmberInpcrdFile(inpcrd_file)
    system = prmtop.createSystem(nonbondedMethod=mmapp.PME, nonbondedCutoff=1*unit.nanometer,
                                 constraints=mmapp.HBonds)
    atm_utils = atm.ATMMetaForceUtils(system)

    number_of_atoms = prmtop.topology.getNumAtoms()

    rcpt_atoms = []
    for at in prmtop.topology.atoms():
        if int(at.residue.id) == rcpt_resid:
            rcpt_atoms.append(at.index)
        
    lig_atoms = []
    for at in prmtop.topology.atoms():
        if int(at.residue.id) == lig_resid:
            lig_atoms.append(at.index)
        
    rcpt_atom_restr = rcpt_atoms
    lig_atom_restr = lig_atoms
     
    kf = 25.0 * unit.kilocalorie_per_mole/unit.angstrom**2 #force constant for Vsite CM-CM restraint
    r0 = 5 * unit.angstrom #radius of Vsite sphere

    #Vsite restraint
    atm_utils.addVsiteRestraintForceCMCM(lig_cm_particles = lig_atom_restr,
                                         rcpt_cm_particles = rcpt_atom_restr,
                                         kfcm = kf,
                                         tolcm = r0,
                                         offset = lig_restr_offset)

    #receptor positional restraints
    fc = 25.0 * unit.kilocalorie_per_mole/unit.angstrom**2
    tol = 0.5 * unit.angstrom
    carbon = re.compile("^C.*")
    posrestr_atoms = []
    for at in prmtop.topology.atoms():
        if int(at.residue.id) == rcpt_resid and carbon.match(at.name) and at.index < 40:
            posrestr_atoms.append(at.index)
    atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)
        
    #create ATM Force
    atmforcegroup = 2
    nonbonded_force_group = 1
    atm_utils.setNonbondedForceGroup(nonbonded_force_group)
    atmvariableforcegroups = [nonbonded_force_group]
    atmforce = atm.ATMMetaForce(lambda1, lambda2,
                                alpha * unit.kilojoules_per_mole, u0/unit.kilojoules_per_mole,
                                w0coeff/unit.kilojoules_per_mole, umsc/unit.kilojoules_per_mole,
                                ubcore/unit.kilojoules_per_mole, acore, direction,
                                atmvariableforcegroups )
    for at in prmtop.topology.atoms():
        atmforce.addParticle(at.index, 0., 0., 0.)
    for i in lig_atoms:
        atmforce.setParticleParameters(i, i, displ[0] * unit.angstrom, displ[1] * unit.angstrom, displ[2] * unit.angstrom)
    atmforce.setForceGroup(atmforcegroup)
    system.addForce(atmforce)
    print("Using ATM Meta Force plugin version = %s" % atm.ATMMETAFORCE_VERSION)

    #set up the integrator
    temperature = 300 * unit.kelvin
    frictionCoeff = 0.5 / unit.picosecond
    MDstepsize = 0.001 * unit.picosecond

    #add barostat, but disabled
    barostat = mm.MonteCarloBarostat(1*unit.bar, temperature)
    barostat.setFrequency(0)#disabled
    system.addForce(barostat)

    integrator = mm.LangevinIntegrator(temperature/unit.kelvin, frictionCoeff/(1/unit.picosecond), MDstepsize/ unit.picosecond)
    integrator.setIntegrationForceGroups({0,atmforcegroup})

    sysplatform = mm.Platform.getPlatformByName(platform)

    properties = {}
    if platform == 'CUDA' or platform == 'OpenCL' or platform == 'HIP':
        properties["Precision"] = "mixed"

    simulation = mmapp.Simulation(prmtop.topology, system, integrator, sysplatform, properties)
    print ("Using platform %s" % simulation.context.getPlatform().getName())
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
    pote = state.getPotentialEnergy()

    print( "LoadState ...")
    simulation.loadState(xml_file)

    #override ATM parameters from checkpoint file
    simulation.context.setParameter(atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(atmforce.Lambda2(), lambda2)
    simulation.context.setParameter(atmforce.Alpha(), alpha * unit.kilojoules_per_mole)
    simulation.context.setParameter(atmforce.U0(), u0 /unit.kilojoules_per_mole)
    simulation.context.setParameter(atmforce.W0(), w0coeff / unit.kilojoules_per_mole)
    simulation.context.setParameter(atmforce.Umax(), umsc / unit.kilojoules_per_mole)
    simulation.context.setParameter(atmforce.Ubcore(), ubcore / unit.kilojoules_per_mole)
    simulation.context.setParameter(atmforce.Acore(), acore)
    simulation.context.setParameter(atmforce.Direction(), acore)

    state = simulation.context.getState(getEnergy = True, groups = {0,atmforcegroup})
    pot_energy = state.getPotentialEnergy()
    pert_energy = atmforce.getPerturbationEnergy(simulation.context)
    print("Potential Energy =", pot_energy)
    print("Binding Energy = ", pert_energy)

    expectedPotEnergy = -116071.0 * unit.kilojoules_per_mole
    expectedPertEnergy = 58.2 * unit.kilojoules_per_mole
    assert np.allclose(expectedPotEnergy.value_in_unit(unit.kilojoules_per_mole), pot_energy.value_in_unit(unit.kilojoules_per_mole), atol = 0.1)
    assert np.allclose(expectedPertEnergy.value_in_unit(unit.kilojoules_per_mole), pert_energy.value_in_unit(unit.kilojoules_per_mole), atol = 0.1)
