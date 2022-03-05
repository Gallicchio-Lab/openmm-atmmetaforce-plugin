#ifndef COMMON_ATMMETAFORCE_KERNELS_H_
#define COMMON_ATMMETAFORCE_KERNELS_H_

#include "ATMMetaForceKernels.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/common/ComputeArray.h"

using namespace OpenMM;

namespace ATMMetaForcePlugin {

/**
 * This kernel is invoked by ATMMetaForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcATMMetaForceKernel : public CalcATMMetaForceKernel {
public:
    CommonCalcATMMetaForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::ComputeContext& cc): CalcATMMetaForceKernel(name, platform), hasInitializedKernel(false), cc(cc) {
    }

    ~CommonCalcATMMetaForceKernel();
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
    
    /**
     * Get the ComputeContext corresponding to the inner Context.
     */
    virtual ComputeContext& getInnerComputeContext(ContextImpl& innerContext) = 0;
    
private:
    class ForceInfo;
    class ReorderListener;
    
    void initkernels(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2);
    
    bool hasInitializedKernel;
    OpenMM::ComputeContext& cc;
    //const OpenMM::System& system;

    //a copy of the displacement vectors stored in the force
    std::vector<mm_float4> displVector;
    
    OpenMM::ComputeArray displ;
    ComputeKernel CopyStateKernel;
    ComputeKernel HybridForceKernel;

    int numParticles;
    double PerturbationEnergy;
};

} // namespace ATMMetaForcePlugin

#endif /*COMMON_ATMMETAFORCE_KERNELS_H_*/
