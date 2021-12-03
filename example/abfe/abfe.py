from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "temoa-g1"

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

rcpt_resid = 1
lig_resid = 2
displ = [ 22.0, 22.0, 22.0 ]

displacement      = [  displ[i] for i in range(3) ] * angstrom
lig_restr_offset  = [    0      for i in range(3) ] * angstrom

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
atm_utils = ATMMetaForceUtils(system)

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
     
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 5 * angstrom #radius of Vsite sphere

#Vsite restraint
atm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                            rcpt_cm_particles = rcpt_atom_restr,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig_restr_offset)

#receptor positional restraints
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
carbon = re.compile("^C.*")
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == rcpt_resid and carbon.match(at.name) and at.index < 40:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)
        
#create ATM Force
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, direction )
for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
for i in lig_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
atmforce.setForceGroup(3)
system.addForce(atmforce)

#set up the integrator
temperature = 300 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

#add barostat, but disabled
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setForceGroup(1)
barostat.setFrequency(0)#disabled
system.addForce(barostat)

integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
integrator.setIntegrationForceGroups({1,3})

platform_name = 'OpenCL'
#platform_name = 'Reference'
platform = Platform.getPlatformByName(platform_name)

properties = {}
properties["Precision"] = "mixed"

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

state = simulation.context.getState(getEnergy = True, groups = {1,3})
pote = state.getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '-equil.xml')

#override ATM parameters from checkpoint file
simulation.context.setParameter(atmforce.Lambda1(), lambda1)
simulation.context.setParameter(atmforce.Lambda2(), lambda2)
simulation.context.setParameter(atmforce.Alpha(), alpha *kilojoules_per_mole)
simulation.context.setParameter(atmforce.U0(), u0 /kilojoules_per_mole)
simulation.context.setParameter(atmforce.W0(), w0coeff /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Umax(), umsc /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Ubcore(), ubcore /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Acore(), acore)
simulation.context.setParameter(atmforce.Direction(), acore)

state = simulation.context.getState(getEnergy = True, groups = {1,3})
print("Potential Energy =", state.getPotentialEnergy())

print("Leg1 production at lambda = %f ..." % lmbd)

stepId = 5000
totalSteps = 50000
loopStep = int(totalSteps/stepId)
simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, potentialEnergy = True, temperature=True))
simulation.reporters.append(DCDReporter(jobname + ".dcd", stepId))

binding_file = jobname + '.out'
f = open(binding_file, 'w')

for i in range(loopStep):
    simulation.step(stepId)
    state = simulation.context.getState(getEnergy = True, groups = {1,3})
    pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
    pert_energy = (atmforce.getPerturbationEnergy(simulation.context)).value_in_unit(kilocalorie_per_mole)
    l1 = simulation.context.getParameter(atmforce.Lambda1())
    l2 = simulation.context.getParameter(atmforce.Lambda2())
    a = simulation.context.getParameter(atmforce.Alpha()) / kilojoules_per_mole
    umid = simulation.context.getParameter(atmforce.U0()) * kilojoules_per_mole
    w0 = simulation.context.getParameter(atmforce.W0()) * kilojoules_per_mole
    print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, l1, l2, a*kilocalorie_per_mole, umid/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, pert_energy), file=f )
    f.flush()

print( "SaveState ...")
simulation.saveState(jobname + '-out.xml')

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
