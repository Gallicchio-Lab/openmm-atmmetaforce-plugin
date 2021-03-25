#ifndef OPENMM_OPENCLATMMETAFORCEKERNELFACTORY_H_
#define OPENMM_OPENCLATMMETAFORCEKERNELFACTORY_H_


#include "openmm/KernelFactory.h"

namespace ATMMetaForcePlugin {

/**
 * This KernelFactory creates kernels for the OpenCL implementation of the ATMMetaForce plugin.
 */

class OpenCLATMMetaForceKernelFactory : public OpenMM::KernelFactory {
public:
    OpenMM::KernelImpl* createKernelImpl(std::string name, const OpenMM::Platform& platform, OpenMM::ContextImpl& context) const;
};

} // namespace ATMMetaForcePlugin

#endif /*OPENMM_OPENCLATMMETAFORCEKERNELFACTORY_H_*/
