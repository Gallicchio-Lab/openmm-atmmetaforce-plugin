#ifndef OPENMM_HIPATMMETAFORCEKERNELFACTORY_H_
#define OPENMM_HIPATMMETAFORCEKERNELFACTORY_H_


#include "openmm/KernelFactory.h"

namespace ATMMetaForcePlugin {

/**
 * This KernelFactory creates kernels for the Hip implementation of the ATMMetaForce plugin.
 */

class HipATMMetaForceKernelFactory : public OpenMM::KernelFactory {
public:
    OpenMM::KernelImpl* createKernelImpl(std::string name, const OpenMM::Platform& platform, OpenMM::ContextImpl& context) const;
};

} // namespace ATMMetaForcePlugin

#endif /*OPENMM_HIPATMMETAFORCEKERNELFACTORY_H_*/
