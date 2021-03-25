#ifndef REFERENCE_ATMMETAFORCE_KERNELS_H_
#define REFERENCE_ATMMETAFORCE_KERNELS_H_

#include "ATMMetaForce.h"
#include "ATMMetaForceKernels.h"
#include "openmm/Platform.h"
#include "openmm/Vec3.h"
#include <vector>

namespace ATMMetaForcePlugin {

/**
 * This kernel is invoked by ATMMetaForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcATMMetaForceKernel : public CalcATMMetaForceKernel {
public:
    ReferenceCalcATMMetaForceKernel(std::string name, const OpenMM::Platform& platform) : CalcATMMetaForceKernel(name, platform) {
    }
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
    double execute(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
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
     * @param innerContext   the context created by the ATM Meta Force for computing displaced energy
     */
    void copyState(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2);


    double getPerturbationEnergy(void) {
      return PerturbationEnergy;
    }
      
 private:
    int numParticles;
    std::vector<int> particles;
    std::vector<OpenMM::Vec3> displ;
    double PerturbationEnergy;

    //softplus parameters
    double lambda1, lambda2, alpha, u0, w0;
    //soft core parameters
    double umax, acore, ubcore;
    
};

} // namespace ATMMetaPlugin

#endif /*REFERENCE_ATMMETAFORCE_KERNELS_H_*/
