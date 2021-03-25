
#include "OpenCLATMMetaForceKernels.h"
#include "OpenCLATMMetaForceKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include <cmath>

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

//NOT SURE ABOUT THIS
class OpenCLATMMetaForceInfo : public OpenCLForceInfo {
public:
    OpenCLATMMetaForceInfo(const ATMMetaForce& force) : OpenCLForceInfo(0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumParticles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        particles.push_back(index);
    }
    bool areGroupsIdentical(int group1, int group2) {
      return (group1 == group2);
    }
private:
    const ATMMetaForce& force;
};

OpenCLCalcATMMetaForceKernel::~OpenCLCalcATMMetaForceKernel() {
}

void OpenCLCalcATMMetaForceKernel::initialize(const System& system, const ATMMetaForce& force) {

  numParticles = force.getNumParticles();
  if (numParticles == 0)
    return;
  vector<mm_float4> displVector(cl.getPaddedNumAtoms());
  for (int i = 0; i < cl.getPaddedNumAtoms(); i++){
    displVector[i].x = 0;
    displVector[i].y = 0;
    displVector[i].z = 0;
    displVector[i].w = 0;
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
  displ = OpenCLArray::create<mm_float4>(cl, cl.getPaddedNumAtoms(), "displ");
  displ->upload(displVector);

  //soft core parameters
  umax = force.getUmax();
  acore = force.getAcore();
  ubcore = force.getUbcore();

  //softplus parameters
  lambda1 = force.getLambda1();
  lambda2 = force.getLambda2();
  alpha = force.getAlpha();
  u0 = force.getU0();
  w0 = force.getW0();
  
  cl.addForce(new OpenCLATMMetaForceInfo(force));
}

void OpenCLCalcATMMetaForceKernel::initkernels(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2){
  if(! hasInitializedKernel) {
    OpenCLContext& cl1 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext1.getPlatformData())->contexts[0];
    OpenCLContext& cl2 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext2.getPlatformData())->contexts[0];

    cl::Program program = cl.createProgram(OpenCLATMMetaForceKernelSources::atmmetaforce, "");
    CopyStateKernel = cl::Kernel(program, "CopyState");
    CopyStateKernel.setArg<cl_int>(0, numParticles);
    CopyStateKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(2, cl1.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(3, cl2.getPosq().getDeviceBuffer());
    CopyStateKernel.setArg<cl::Buffer>(4, displ->getDeviceBuffer());

    HybridForceKernel = cl::Kernel(program, "HybridForce");
    HybridForceKernel.setArg<cl_int>(0, numParticles);
    HybridForceKernel.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
    HybridForceKernel.setArg<cl::Buffer>(2, cl1.getForce().getDeviceBuffer());
    HybridForceKernel.setArg<cl::Buffer>(3, cl2.getForce().getDeviceBuffer());
    //there is a fifth argument (sp) which is added in execute()

    hasInitializedKernel = true;
  }
}

double OpenCLCalcATMMetaForceKernel::execute(OpenMM::ContextImpl& context,
					     OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
					     double State1Energy, double State2Energy,
					     bool includeForces, bool includeEnergy ) {
  initkernels(context, innerContext1, innerContext2);

  //soft-core perturbation energy
  double fp;
  PerturbationEnergy = SoftCoreF(State2Energy - State1Energy, umax, acore, ubcore, fp);

  //softplus function
  double ebias = 0.0;
  double ee = 1.0 + exp(-alpha*(PerturbationEnergy  - u0));
  if(alpha > 0){
    ebias = ((lambda2 - lambda1)/alpha) * log(ee);
  }
  ebias += lambda2 * PerturbationEnergy  + w0;
  double bfp = (lambda2 - lambda1)/ee + lambda1;

  //alchemical potential energy
  double energy = State1Energy + ebias;

  //hybridize forces and add them to the system's forces
  float sp = bfp*fp;
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
  //if (force.getNumParticles() != particles.size())
  //      throw OpenMMException("copyParametersToContext: The number of ATMMetaForce particles has changed");
  OpenCLContext& cl = *reinterpret_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
  vector<mm_float4> displVector(cl.getPaddedNumAtoms());
  for (int i = 0; i < numParticles; i++){
    int particle;
    double dx, dy, dz;
    force.getParticleParameters(i, particle, dx, dy, dz);
    displVector[i].x = dx;
    displVector[i].y = dy;
    displVector[i].z = dz;
    displVector[i].w = 0;
  }
  
  //soft core parameters
  umax = force.getUmax();
  acore = force.getAcore();
  ubcore = force.getUbcore();

  //softplus parameters
  lambda1 = force.getLambda1();
  lambda2 = force.getLambda2();
  alpha = force.getAlpha();
  u0 = force.getU0();
  w0 = force.getW0();
}
