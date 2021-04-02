from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from desmonddmsfile75 import *
from atmmetaforce import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

binding_file = 'temoa-g1-g4.out'
f = open(binding_file, 'w')

temperature = 300.0 * kelvin
lmbd = 0.50
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole

rcpt_resid = 1
lig1_resid = 2
lig2_resid = 3
displ = [22.0, 22.0, 22.0]

displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

refatoms_lig1 = [8, 6, 4]
refatoms_lig2 = [3, 5, 1]


umsc =  100.0 * kilocalorie_per_mole
ubcore = 50.0 * kilocalorie_per_mole
acore = 0.062500

print("temperature = ", temperature)
print("lambda = ", lmbd)
print("lambda1 = ", lambda1)
print("lambda2 = ", lambda2)
print("alpha = ", alpha)
print("u0 = ", u0)
print("w0coeff = ", w0coeff)
print("soft core method = ", 'RationalSoftCoreMethod')
print("umax = ", umsc)
print("acore = ", acore)
print("ubcore = ", ubcore)
print("displacement = ", displacement)

file_input  = 'temoa-g1-g4.dms'
file_output = 'temoa-g1-g4-out.dms'

shutil.copyfile(file_input, file_output)

testDes = DesmondDMSFile(file_output) 
system = testDes.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*nanometer)

number_of_atoms = testDes.topology.getNumAtoms()

lig1_atoms = []
for at in testDes.topology.atoms():
    if int(at.residue.id) == lig1_resid:
        lig1_atoms.append(int(at.id)-1)
        
lig2_atoms = []
for at in testDes.topology.atoms():
    if int(at.residue.id) == lig2_resid:
        lig2_atoms.append(int(at.id)-1)

lig1_atom_restr = lig1_atoms
lig2_atom_restr = lig2_atoms

rcpt_atom_restr = []
for at in testDes.topology.atoms():
    if int(at.residue.id) == rcpt_resid:
        rcpt_atom_restr.append(int(at.id)-1)
     
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 4.5 * angstrom #radius of Vsite sphere

#these can be 'None" if not using orientational restraints
lig_ref_atoms = None # the 3 atoms of the ligand that define the coordinate system of the ligand
rcpt_ref_atoms = None # the 3 atoms of the receptor that define the coordinate system of the receptor
angle_center = None * degrees
kfangle = None * kilocalorie_per_mole/degrees**2
angletol = None * degrees
dihedral1center = None * degrees
kfdihedral1 = None * kilocalorie_per_mole/degrees**2
dihedral1tol = None * degrees
dihedral2center = None * degrees
kfdihedral2 = None * kilocalorie_per_mole/degrees**2
dihedral2tol = None * degrees

sdm_utils = ATMMetaForceUtils(system)
sdm_utils.addRestraintForce(lig_cm_particles = lig1_atom_restr,
                            rcpt_cm_particles = rcpt_atom_restr,
                            kfcm = kf,
                            tolcm = r0,
                            lig_ref_particles = lig_ref_atoms,
                            rcpt_ref_particles = rcpt_ref_atoms,
                            angle_center = angle_center,
                            kfangle = kfangle,
                            angletol = angletol,
                            dihedral1center = dihedral1center,
                            kfdihedral1 = kfdihedral1,
                            dihedral1tol = dihedral1tol,
                            dihedral2center = dihedral2center,
                            kfdihedral2 = kfdihedral2,
                            dihedral2tol = dihedral2tol,
                            offset = lig1_restr_offset)

sdm_utils.addRestraintForce(lig_cm_particles = lig2_atom_restr,
                            rcpt_cm_particles = rcpt_atom_restr,
                            kfcm = kf,
                            tolcm = r0,
                            lig_ref_particles = lig_ref_atoms,
                            rcpt_ref_particles = rcpt_ref_atoms,
                            angle_center = angle_center,
                            kfangle = kfangle,
                            angletol = angletol,
                            dihedral1center = dihedral1center,
                            kfdihedral1 = kfdihedral1,
                            dihedral1tol = dihedral1tol,
                            dihedral2center = dihedral2center,
                            kfdihedral2 = kfdihedral2,
                            dihedral2tol = dihedral2tol,
                            offset = lig2_restr_offset)

lig1_ref_atoms  = [ refatoms_lig1[i]+lig1_atoms[0] for i in range(3)]
lig2_ref_atoms  = [ refatoms_lig2[i]+lig2_atoms[0] for i in range(3)]
sdm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                            ligb_ref_particles = lig2_ref_atoms,
                            kfdispl = 25.0 * kilocalorie_per_mole/angstrom**2,
                            ktheta =  50.0 * kilocalorie_per_mole,
                            kpsi =  50.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)


#platform_name = 'Reference'
platform_name = 'OpenCL'
platform = Platform.getPlatformByName(platform_name)

properties = {}

if platform_name =='OpenCL':
    #expected "platformid:deviceid" or empty
    #device = "@pn@"
    device = "0:0"
    m = re.match("(\d+):(\d+)", device)
    if m:
        platformid = m.group(1)
        deviceid = m.group(2)
        properties["OpenCLPlatformIndex"] = platformid
        properties["DeviceIndex"] = deviceid
        print("Using platform id: %s, device id: %s" % ( platformid ,  deviceid) )

#create ATM Force
atmforce = ATMMetaForce()

for at in testDes.topology.atoms():
    atmforce.addParticle(int(at.id)-1, 0., 0., 0.)

for i in lig1_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
    
for i in lig2_atoms:
    atmforce.setParticleParameters(i, i, -displ[0] * angstrom, -displ[1] * angstrom, -displ[2] * angstrom)

atmforce.setLambda1(lambda1);
atmforce.setLambda2(lambda2);
atmforce.setUmax(umsc/kilojoules_per_mole);
atmforce.setUbcore(ubcore/kilojoules_per_mole);
atmforce.setAcore(acore);
    
atmforce.setForceGroup(3)

system.addForce(atmforce)



#add Langevin integrator   
frictionCoeff = 0.1 / picosecond
MDstepsize = 0.001 * picosecond

integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond )
integrator.setIntegrationForceGroups({1,3})

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)

state = simulation.context.getState(getEnergy = True, groups = {1})
print("Group 1: ", state.getPotentialEnergy())

state = simulation.context.getState(getEnergy = True, groups = {3})
print("Group 3: ", state.getPotentialEnergy())

totalSteps = 10000
nprnt = 1000
ntrj = 1000
simulation.reporters.append(StateDataReporter(stdout, nprnt, step=True, temperature=True))
simulation.reporters.append(DCDReporter("temoa-g1-g4.dcd", ntrj))

loops = int(totalSteps/nprnt)
start=datetime.now()
step = 0

for i in range(loops):
    simulation.step(nprnt)
    state = simulation.context.getState(getEnergy = True, groups = {1,3})
    pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
    pert_energy = (atmforce.getPerturbationEnergy(simulation.context)).value_in_unit(kilocalorie_per_mole)
    print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, lambda1, lambda2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0coeff/kilocalorie_per_mole, pot_energy, pert_energy), file=f )
    f.flush()
    step += nprnt
end=datetime.now()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

f.close()

elapsed=end - start
print("MD time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
