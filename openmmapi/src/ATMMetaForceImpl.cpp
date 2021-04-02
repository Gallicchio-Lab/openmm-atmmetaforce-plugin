
#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "internal/ATMMetaForceImpl.h"
#include "ATMMetaForceKernels.h"

#include "openmm/NonbondedForce.h"
#include "openmm/kernels.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/Vec3.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include <iostream>

using namespace ATMMetaForcePlugin;
using namespace OpenMM;
using namespace std;

ATMMetaForceImpl::ATMMetaForceImpl(const ATMMetaForce& owner) : owner(owner), innerIntegrator1(1.0), innerIntegrator2(1.0), hasInitializedInnerContexts(false){
}

ATMMetaForceImpl::~ATMMetaForceImpl() {
}

static void copysystem(const OpenMM::System& system, OpenMM::System& innerSystem){
  
  //copy particles
  for (int i = 0; i < system.getNumParticles(); i++) 
    innerSystem.addParticle(system.getParticleMass(i));

  //copy constraints
  for (int i = 0; i < system.getNumConstraints(); i++) {
    int particle1, particle2;
    double distance;
    system.getConstraintParameters(i, particle1, particle2, distance);
    innerSystem.addConstraint(particle1, particle2, distance);
  }

  Vec3 a, b, c;
  system.getDefaultPeriodicBoxVectors(a, b, c);
  innerSystem.setDefaultPeriodicBoxVectors(a, b, c);
  
  //add non-bonded force
  int numforces = system.getNumForces();
  for (int i=0; i<numforces; i++){
    const Force &force = system.getForce(i);
    int group = force.getForceGroup();
    if(group == 2){//non-bonded
      Force* newnbforce = XmlSerializer::clone<Force>(force);
      newnbforce->setForceGroup(group);
      NonbondedForce* nonbonded = dynamic_cast<NonbondedForce*>(newnbforce);
      if (nonbonded != NULL)
	nonbonded->setReciprocalSpaceForceGroup(-1);
      innerSystem.addForce(newnbforce);
    }
  }

}
  
void ATMMetaForceImpl::initialize(ContextImpl& context) {
  const OpenMM::System& system = context.getSystem();

  int numforces = system.getNumForces();
  for (int i=0; i<numforces; i++){
    const Force &force = system.getForce(i);
    int group = force.getForceGroup();
    if (!(group == 1 || group == 2 || group == 3)){
      throw OpenMMException("The ATM Meta Force requires all forces to be either in force group 1 (bonded), group 2 (nonbonded) or group 3 (ATM Meta Force)");
    }
  }

  copysystem(system, innerSystem1);
  copysystem(system, innerSystem2);

  kernel = context.getPlatform().createKernel(CalcATMMetaForceKernel::Name(), context);
  kernel.getAs<CalcATMMetaForceKernel>().initialize(context.getSystem(), owner);
}

double ATMMetaForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
  if(! hasInitializedInnerContexts) {
    hasInitializedInnerContexts = true;
    
    innerContext1 = context.createLinkedContext(innerSystem1, innerIntegrator1);
    innerContext2 = context.createLinkedContext(innerSystem2, innerIntegrator2);

    vector<Vec3> pos;
    context.getPositions(pos);
    innerContext1->setPositions(pos);
    innerContext2->setPositions(pos);
  }

  
  if ((groups&(1<<owner.getForceGroup())) == 0) return 0.0;
  
  bool do_energy = true;
  ContextImpl& innercontextimpl1 = getContextImpl(*innerContext1);
  ContextImpl& innercontextimpl2 = getContextImpl(*innerContext2);

  //copies the coordinates etc. from the context to the inner contexts 
  kernel.getAs<CalcATMMetaForceKernel>().copyState(context, innercontextimpl1, innercontextimpl2);
 
  //evaluate energy and force for original system, 4 = evaluate force group 2 (non-bonded)
  double State1Energy = innercontextimpl1.calcForcesAndEnergy(true, do_energy, 4);

  //evaluate energy and force for the displaced system, 4 = evaluate force group 2 (non-bonded)
  double State2Energy = innercontextimpl2.calcForcesAndEnergy(true, do_energy, 4);

  //evalaute the alchemical energy
  double energy = kernel.getAs<CalcATMMetaForceKernel>().execute(context,
								 innercontextimpl1, innercontextimpl2,
								 State1Energy, State2Energy,
								 includeForces, includeEnergy);

  //retrieve the perturbation energy
  PerturbationEnergy = kernel.getAs<CalcATMMetaForceKernel>().getPerturbationEnergy();

  return (includeEnergy ? energy : 0.0);
}

std::vector<std::string> ATMMetaForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcATMMetaForceKernel::Name());
    return names;
}

void ATMMetaForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcATMMetaForceKernel>().copyParametersToContext(context, owner);
}
