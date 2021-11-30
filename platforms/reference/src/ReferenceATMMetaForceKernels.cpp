
#include "ReferenceATMMetaForceKernels.h"
#include "ATMMetaForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Vec3.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <iostream>
#include <cmath>

using namespace ATMMetaForcePlugin;
using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

static double SoftCoreF(double u, double umax, double a, double ub, double& fp){
  if(u <= ub){
    fp = 1.;
    return u;
  }
  double gu = (u-ub)/(a*(umax-ub)); //this is y/alpha
  double zeta = 1. + 2.*gu*(gu + 1.) ;
  double zetap = pow( zeta , a );
  double s = 4.*(2.*gu + 1.)/zeta;
  fp = s*zetap/pow(1.+zetap,2);
  return (umax-ub)*(zetap - 1.)/(zetap + 1.) + ub;
}



void ReferenceCalcATMMetaForceKernel::initialize(const System& system, const ATMMetaForce& force) {

  int numParticles = force.getNumParticles();
  particles.resize(numParticles);
  
  //displacement map
  displ.resize(numParticles);
  for (int i = 0; i < numParticles; i++){
    double dx, dy, dz;
    force.getParticleParameters(i, particles[i], dx, dy, dz);
    displ[i] = OpenMM::Vec3(dx,dy,dz);
  }
  
}

double ReferenceCalcATMMetaForceKernel::execute(ContextImpl& context, ContextImpl& innerContext1, ContextImpl& innerContext2,
						double State1Energy, double State2Energy,
						bool includeForces, bool includeEnergy) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);
    vector<Vec3>& force1 = extractForces(innerContext1);
    vector<Vec3>& force2 = extractForces(innerContext2);
    int numParticles = particles.size();

    //softplus parameters
    double lambda1 = context.getParameter(ATMMetaForce::Lambda1());
    double lambda2 = context.getParameter(ATMMetaForce::Lambda2());
    double alpha   = context.getParameter(ATMMetaForce::Alpha());
    double u0      = context.getParameter(ATMMetaForce::U0());
    double w0      = context.getParameter(ATMMetaForce::W0());

    //softcore parameters
    double umax = context.getParameter(ATMMetaForce::Umax());
    double ubcore = context.getParameter(ATMMetaForce::Ubcore());
    double acore = context.getParameter(ATMMetaForce::Acore());

    //alchemical direction
    // 1 = from RA  (reference) to R+A (displaced)
    //-1 = from R+A (displaced) to RA  (reference)
    double alchemical_direction = context.getParameter(ATMMetaForce::Direction());

    //soft-core perturbation energy
    double fp;
    double u  = alchemical_direction > 0 ?  State2Energy - State1Energy : State1Energy - State2Energy;
    double e0 = alchemical_direction > 0 ?  State1Energy                : State2Energy;
    PerturbationEnergy = SoftCoreF(u, umax, acore, ubcore, fp);

    //softplus function
    double ebias = 0.0;
    double ee = 1.0 + exp(-alpha*(PerturbationEnergy  - u0));
    if(alpha > 0){
      ebias = ((lambda2 - lambda1)/alpha) * log(ee);
    }
    ebias += lambda2 * PerturbationEnergy  + w0;
    double bfp = (lambda2 - lambda1)/ee + lambda1;

    //alchemical potential energy
    double energy = e0 + ebias;

    //hybridize forces and add them to the system's forces
    double sp = bfp*fp;
    if(alchemical_direction > 0){
      for(int i=0; i < numParticles; i++){
	force[i] += sp*force2[i] + (1.-sp)*force1[i];
      }
    }else{
      for(int i=0; i < numParticles; i++){
	force[i] += sp*force1[i] + (1.-sp)*force2[i];
      }
    }

    return (includeEnergy ? energy : 0.0);
}


void ReferenceCalcATMMetaForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext1, ContextImpl& innerContext2) {
  vector<Vec3>& pos = extractPositions(context);
  extractPositions(innerContext1) = pos;

  vector<Vec3> pos2(pos);
  for(int i=0; i < pos2.size(); i++){
    pos2[i] += displ[i];
  }
  extractPositions(innerContext2) = pos2;

  Vec3 a, b, c;
  context.getPeriodicBoxVectors(a, b, c);
  innerContext1.setPeriodicBoxVectors(a, b, c);
  innerContext2.setPeriodicBoxVectors(a, b, c);
  
  innerContext1.setTime(context.getTime());
  innerContext2.setTime(context.getTime());
  
  map<string, double> innerParameters;

  innerParameters = innerContext1.getParameters();
  for (auto& param : innerParameters)
    innerContext1.setParameter(param.first, context.getParameter(param.first));
  
  innerParameters = innerContext2.getParameters();
  for (auto& param : innerParameters)
    innerContext2.setParameter(param.first, context.getParameter(param.first));
  
}

void ReferenceCalcATMMetaForceKernel::copyParametersToContext(ContextImpl& context, const ATMMetaForce& force) {
    if (force.getNumParticles() != particles.size())
        throw OpenMMException("copyParametersToContext: The number of ATMMetaForce particles has changed");
    for (int i = 0; i < force.getNumParticles(); i++) {
      int p;
      double dx, dy, dz;
      force.getParticleParameters(i, p, dx, dy, dz);
      if (p != particles[i])
	throw OpenMMException("ReferenceCalcATMMetaForceKernel::copyParametersToContext: A particle index has changed");
      displ.push_back(OpenMM::Vec3(dx,dy,dz));
    }
}
