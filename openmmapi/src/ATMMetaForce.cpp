
#include "ATMMetaForce.h"
#include "internal/ATMMetaForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>

using namespace ATMMetaForcePlugin;
using namespace OpenMM;
using namespace std;

ATMMetaForce::ATMMetaForce() {

  //default values for soft-plus potential
  lambda1 = lambda2 = 1.0;
  alpha = 1.0;
  u0 = 0.0;

  //default values for the soft-core function
  umax = 200.0;
  acore = 0.25;
  ubcore = 0.0;
  
}

int ATMMetaForce::addParticle(int particle, double dx, double dy, double dz) {
    particles.push_back(ParticleInfo(particle, dx, dy, dz));
    return particles.size()-1;
}

void ATMMetaForce::getParticleParameters(int index, int& particle, double& dx, double &dy, double &dz) const {
    ASSERT_VALID_INDEX(index, particles);
    particle = particles[index].particle;
    dx = particles[index].dx;
    dy = particles[index].dy;
    dz = particles[index].dz;

}

void ATMMetaForce::setParticleParameters(int index, int particle, double dx, double dy, double dz){
    ASSERT_VALID_INDEX(index, particles);
    particles[index].particle = particle;
    particles[index].dx = dx;
    particles[index].dy = dy;
    particles[index].dz = dz;
}

ForceImpl* ATMMetaForce::createImpl() const {
    return new ATMMetaForceImpl(*this);
}

void ATMMetaForce::updateParametersInContext(OpenMM::Context& context) {
  dynamic_cast<ATMMetaForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

double ATMMetaForce::getPerturbationEnergy(const OpenMM::Context& context)  const {
  return dynamic_cast<const ATMMetaForceImpl&>(getImplInContext(context)).getPerturbationEnergy(); 
}

