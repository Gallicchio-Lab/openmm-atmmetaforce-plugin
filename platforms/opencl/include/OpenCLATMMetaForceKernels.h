#ifndef OPENCL_ATMMETAFORCE_KERNELS_H_
#define OPENCL_ATMMETAFORCE_KERNELS_H_

#include "ATMMetaForceKernels.h"
#include "openmm/opencl/OpenCLContext.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLArray.h"

using namespace OpenMM;

namespace ATMMetaForcePlugin {

/**
 * This kernel is invoked by ATMMetaForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcATMMetaForceKernel : public CalcATMMetaForceKernel {
public:
    OpenCLCalcATMMetaForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::OpenCLContext& cl, const OpenMM::System& system) :
            CalcATMMetaForceKernel(name, platform), hasInitializedKernel(false), cl(cl), system(system), displ(NULL) {
    }

    ~OpenCLCalcATMMetaForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ATMMetaForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const ATMMetaForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context,
		   OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
		   double State1Energy, double State2Energy,
		   bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ATMMetaForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const ATMMetaForce& force);

    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     */
    void copyState(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2);


    double getPerturbationEnergy(void) {
      return PerturbationEnergy;
    }
    
private:
    class ForceInfo;
    class ReorderListener;
    
    void initkernels(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2);
    
    bool hasInitializedKernel;
    OpenMM::OpenCLContext& cl;
    const OpenMM::System& system;

    std::vector<mm_float4> displVector;
    
    OpenMM::OpenCLArray* displ;
    cl::Kernel CopyStateKernel;
    cl::Kernel HybridForceKernel;

    int numParticles;
    double PerturbationEnergy;

    //softplus parameters
    double lambda1, lambda2, alpha, u0, w0;
    //soft core parameters
    double umax, acore, ubcore;
    
};

} // namespace ATMMetaForcePlugin

#endif /*OPENCL_ATMMETAFORCE_KERNELS_H_*/
