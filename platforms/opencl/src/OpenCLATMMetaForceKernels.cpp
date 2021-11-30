#include "ATMMetaForce.h"
#include "OpenCLATMMetaForceKernels.h"
#include "OpenCLATMMetaForceKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include <cmath>

//DEBUG
#include <iostream>

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

//NOT SURE ABOUT THIS, it's not added to the Context?
class OpenCLCalcATMMetaForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(ComputeForceInfo& force) : OpenCLForceInfo(0), force(force) {
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


class OpenCLCalcATMMetaForceKernel::ReorderListener : public OpenCLContext::ReorderListener {
public:
  ReorderListener(OpenCLContext&  cl, vector<mm_float4>& displVector, OpenCLArray* displ) :
    cl(cl), displVector(displVector), displ(displ)  {
    }
    void execute() {
        const vector<int>& id = cl.getAtomIndex();
	vector<mm_float4> newDisplVectorContext(cl.getPaddedNumAtoms());
	for (int i = 0; i < cl.getNumAtoms(); i++){
	  newDisplVectorContext[i] = displVector[id[i]];
	}
	displ->upload(newDisplVectorContext);
    }
private:
    OpenCLContext& cl;
    OpenCLArray* displ;
    std::vector<mm_float4> displVector;
  
  
};


OpenCLCalcATMMetaForceKernel::~OpenCLCalcATMMetaForceKernel() {
}

void OpenCLCalcATMMetaForceKernel::initialize(const System& system, const ATMMetaForce& force) {

  numParticles = force.getNumParticles();
  if (numParticles == 0)
    return;
  displVector.resize(cl.getPaddedNumAtoms());
  vector<mm_float4> displVectorContext(cl.getPaddedNumAtoms());
  for (int i = 0; i < cl.getPaddedNumAtoms(); i++){
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
  const vector<int>& id = cl.getAtomIndex();
  for (int i = 0; i < numParticles; i++){
    displVectorContext[i] = displVector[id[i]];
  }
  displ = OpenCLArray::create<mm_float4>(cl, cl.getPaddedNumAtoms(), "displ");
  displ->upload(displVectorContext);

  cl.addForce(new OpenCLForceInfo(1));
}

void OpenCLCalcATMMetaForceKernel::initkernels(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2){
  if(! hasInitializedKernel) {
    hasInitializedKernel = true;

    //inner contexts
    OpenCLContext& cl1 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext1.getPlatformData())->contexts[0];
    OpenCLContext& cl2 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext2.getPlatformData())->contexts[0];

    //initialize the listener, this reorders the displacement vectors
    ReorderListener* listener = new ReorderListener(cl, displVector, displ );
    cl.addReorderListener(listener);
    listener->execute();

    //create CopyState kernel
    cl::Program program = cl.createProgram(OpenCLATMMetaForceKernelSources::atmmetaforce, "");
    CopyStateKernel = cl::Kernel(program, "CopyState");
    CopyStateKernel.setArg<cl_int>(0, numParticles);
    int nargs = 0;
    CopyStateKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(2, cl1.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(3, cl2.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(4, displ->getDeviceBuffer());
    if(cl.getUseMixedPrecision()){
      CopyStateKernel.setArg<cl::Buffer>(5, cl.getPosqCorrection().getDeviceBuffer());
      CopyStateKernel.setArg<cl::Buffer>(6, cl1.getPosqCorrection().getDeviceBuffer());
      CopyStateKernel.setArg<cl::Buffer>(7, cl2.getPosqCorrection().getDeviceBuffer());
    }

    //create the HybridForce kernel
    HybridForceKernel = cl::Kernel(program, "HybridForce");
    HybridForceKernel.setArg<cl_int>(0, numParticles);
    HybridForceKernel.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
    HybridForceKernel.setArg<cl::Buffer>(2, cl1.getForce().getDeviceBuffer());
    HybridForceKernel.setArg<cl::Buffer>(3, cl2.getForce().getDeviceBuffer());
    //there is a 4th argument (sp) which is added in execute()

    cl1.addForce(new OpenCLForceInfo(1));
    cl2.addForce(new OpenCLForceInfo(1));
  }
}

double OpenCLCalcATMMetaForceKernel::execute(OpenMM::ContextImpl& context,
					     OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
					     double State1Energy, double State2Energy,
					     bool includeForces, bool includeEnergy ) {
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
  HybridForceKernel.setArg<cl_float>(4, sp);
  cl.executeKernel(HybridForceKernel, numParticles);
  
  return (includeEnergy ? energy : 0.0);
}

void OpenCLCalcATMMetaForceKernel::copyState(OpenMM::ContextImpl& context,
					     OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2) {
  initkernels(context, innerContext1, innerContext2);
  cl.executeKernel(CopyStateKernel, numParticles);
  
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


void OpenCLCalcATMMetaForceKernel::copyParametersToContext(ContextImpl& context, const ATMMetaForce& force) {
  OpenCLContext& cl = *reinterpret_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
  for (int i = 0; i < numParticles; i++){
    int particle;
    double dx, dy, dz;
    force.getParticleParameters(i, particle, dx, dy, dz);
    displVector[i].x = dx;
    displVector[i].y = dy;
    displVector[i].z = dz;
    displVector[i].w = 0;
  }
  const vector<int>& id = cl.getAtomIndex();
  vector<mm_float4> displVectorContext(displVector);
  for (int i = 0; i < cl.getPaddedNumAtoms(); i++){
    displVectorContext[i].x = 0;
    displVectorContext[i].y = 0;
    displVectorContext[i].z = 0;
    displVectorContext[i].w = 0;
  }
  for (int i = 0; i < numParticles; i++){
    displVectorContext[i] = displVector[id[i]];
  }
  displ->upload(displVectorContext);
}
