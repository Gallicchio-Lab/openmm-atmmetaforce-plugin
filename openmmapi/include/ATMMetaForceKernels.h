#ifndef ATMMETAFORCE_KERNELS_H_
#define ATMMETAFORCE_KERNELS_H_

#include "ATMMetaForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace ATMMetaForcePlugin {

/**
 * This kernel is invoked by ATMMetaForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcATMMetaForceKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "CalcATMMetaForce";
    }
    CalcATMMetaForceKernel(std::string name, const OpenMM::Platform& platform) : OpenMM::KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ATMMetaForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const ATMMetaForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2,
			   double State1Energy, double State2Energy,
			   bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ATMMetaForce to copy the parameters from
     */
    virtual void copyParametersToContext(OpenMM::ContextImpl& context, const ATMMetaForce& force) = 0;


    /**
     * Copy state information to the inner context.
     *
     * @param context        the context in which to execute this kernel
     * @param innerContext   the context created by the ATM Meta Force for computing displaced energy
     */
    virtual void copyState(OpenMM::ContextImpl& context, OpenMM::ContextImpl& innerContext1, OpenMM::ContextImpl& innerContext2 ) = 0;

    virtual double getPerturbationEnergy(void) = 0;

};



} // namespace ATMMetaPlugin

#endif /*ATMMETAFORCE_KERNELS_H_*/
