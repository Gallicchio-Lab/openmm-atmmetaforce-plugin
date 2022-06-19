#include "ATMMetaForce.h"
#include "CommonATMMetaForceKernels.h"
#include "CommonATMMetaForceKernelSources.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ContextSelector.h" // requires OpenMM 7.7+
#include <cmath>

//DEBUG
//#include <iostream>

using namespace ATMMetaForcePlugin;
using namespace OpenMM;
using namespace std;

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

class CommonCalcATMMetaForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(ComputeForceInfo& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return force.areParticlesIdentical(particle1, particle2);
    }
    int getNumParticleGroups() {
        return force.getNumParticleGroups();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        force.getParticlesInGroup(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return force.areGroupsIdentical(group1, group2);
    }
private:
    ComputeForceInfo& force;
};


class CommonCalcATMMetaForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
  ReorderListener(ComputeContext&  cc, vector<mm_float4>& displVector, ArrayInterface& displ) :
    cc(cc), displVector(displVector), displ(displ)  {
    }
    void execute() {
        const vector<int>& id = cc.getAtomIndex();
	vector<mm_float4> newDisplVectorContext(cc.getPaddedNumAtoms());
	for (int i = 0; i < cc.getNumAtoms(); i++){
	  newDisplVectorContext[i] = displVector[id[i]];
	}
	displ.upload(newDisplVectorContext);
    }
private:
    ComputeContext& cc;
    ArrayInterface& displ;
    std::vector<mm_float4> displVector;
  
  
};


CommonCalcATMMetaForceKernel::~CommonCalcATMMetaForceKernel() {
}

void CommonCalcATMMetaForceKernel::initialize(const System& system, const ATMMetaForce& force) {
  ContextSelector selector(cc); //requires OpenMM 7.7+
  numParticles = force.getNumParticles();
  if (numParticles == 0)
    return;
  displVector.resize(cc.getPaddedNumAtoms());
  vector<mm_float4> displVectorContext(cc.getPaddedNumAtoms());
  for (int i = 0; i < cc.getPaddedNumAtoms(); i++){
    displVector[i].x = displVectorContext[i].x = 0;
    displVector[i].y = displVectorContext[i].y = 0;
    displVector[i].z = displVectorContext[i].z = 0;
    displVector[i].w = displVectorContext[i].w = 0;
  }
  for (int i = 0; i < numParticles; i++){
    int particle;
    double dx, dy, dz;
    force.getParticleParameters(i, particle, dx, dy, dz);
    displVector[i].x = dx;
    displVector[i].y = dy;
    displVector[i].z = dz;
    displVector[i].w = 0;
  }
  const vector<int>& id = cc.getAtomIndex();
  for (int i = 0; i < numParticles; i++){
    displVectorContext[i] = displVector[id[i]];
  }
  displ = ComputeArray();
  displ.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "displ");
  displ.upload(displVectorContext);

  cc.addForce(new ComputeForceInfo());
}

void CommonCalcATMMetaForceKernel::initkernels(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2){
  if(! hasInitializedKernel) {
    hasInitializedKernel = true;

    //inner contexts
    ComputeContext& cc1 = getInnerComputeContext(innerContext1);
    ComputeContext& cc2 = getInnerComputeContext(innerContext2);

    //initialize the listener, this reorders the displacement vectors
    ReorderListener* listener = new ReorderListener(cc, displVector, displ );
    cc.addReorderListener(listener);
    listener->execute();

    //create CopyState kernel
    ComputeProgram program = cc.compileProgram(CommonATMMetaForceKernelSources::atmmetaforce);
    CopyStateKernel = program->createKernel("CopyState");
    CopyStateKernel->addArg(numParticles);
    CopyStateKernel->addArg( cc.getPosq());
    CopyStateKernel->addArg(cc1.getPosq());
    CopyStateKernel->addArg(cc2.getPosq());
    CopyStateKernel->addArg(displ);
    if(cc.getUseMixedPrecision()){
      CopyStateKernel->addArg( cc.getPosqCorrection());
      CopyStateKernel->addArg(cc1.getPosqCorrection());
      CopyStateKernel->addArg(cc2.getPosqCorrection());
    }

    //create the HybridForce kernel
    float sp = 0;
    HybridForceKernel = program->createKernel("HybridForce");
    HybridForceKernel->addArg(numParticles);
    HybridForceKernel->addArg(cc.getPaddedNumAtoms());
    HybridForceKernel->addArg( cc.getLongForceBuffer());
    HybridForceKernel->addArg(cc1.getLongForceBuffer());
    HybridForceKernel->addArg(cc2.getLongForceBuffer());
    HybridForceKernel->addArg(sp);//argument 5 (sp) is set in execute()

    cc1.addForce(new ComputeForceInfo());
    cc2.addForce(new ComputeForceInfo());
    
  }
}

double CommonCalcATMMetaForceKernel::execute(OpenMM::ContextImpl& context,
					     OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
					     double State1Energy, double State2Energy,
					     bool includeForces, bool includeEnergy ) {
  ContextSelector selector(cc); // requires OpenMM 7.7+
  //cc.setAsCurrent(); //OpenMM 7.6

  initkernels(context, innerContext1, innerContext2);

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
  float sp = alchemical_direction > 0 ? bfp*fp : 1. - bfp*fp;
  HybridForceKernel->setArg(5, sp);
  HybridForceKernel->execute(numParticles);
  
  return (includeEnergy ? energy : 0.0);
}

void CommonCalcATMMetaForceKernel::copyState(OpenMM::ContextImpl& context,
					     OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2) {
  ContextSelector selector(cc); // requires OpenMM 7.7

  initkernels(context, innerContext1, innerContext2);

  CopyStateKernel->execute(numParticles);
  
  Vec3 a, b, c;
  context.getPeriodicBoxVectors(a, b, c);
  innerContext1.setPeriodicBoxVectors(a, b, c);
  innerContext1.setTime(context.getTime());
  innerContext2.setPeriodicBoxVectors(a, b, c);
  innerContext2.setTime(context.getTime());
  map<string, double> innerParameters1 = innerContext1.getParameters();
  for (auto& param : innerParameters1)
    innerContext1.setParameter(param.first, context.getParameter(param.first));
  map<string, double> innerParameters2 = innerContext2.getParameters();
  for (auto& param : innerParameters2)
    innerContext2.setParameter(param.first, context.getParameter(param.first));
}


void CommonCalcATMMetaForceKernel::copyParametersToContext(ContextImpl& context, const ATMMetaForce& force) {
  for (int i = 0; i < numParticles; i++){
    int particle;
    double dx, dy, dz;
    force.getParticleParameters(i, particle, dx, dy, dz);
    displVector[i].x = dx;
    displVector[i].y = dy;
    displVector[i].z = dz;
    displVector[i].w = 0;
  }
  const vector<int>& id = cc.getAtomIndex();
  vector<mm_float4> displVectorContext(displVector);
  for (int i = 0; i < cc.getPaddedNumAtoms(); i++){
    displVectorContext[i].x = 0;
    displVectorContext[i].y = 0;
    displVectorContext[i].z = 0;
    displVectorContext[i].w = 0;
  }
  for (int i = 0; i < numParticles; i++){
    displVectorContext[i] = displVector[id[i]];
  }
  displ.upload(displVectorContext);
}
