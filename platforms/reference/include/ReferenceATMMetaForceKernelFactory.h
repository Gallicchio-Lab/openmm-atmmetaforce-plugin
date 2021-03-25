#ifndef OPENMM_REFERENCEATMMETAFORCEKERNELFACTORY_H_
#define OPENMM_REFERENCEATMMETAFORCEKERNELFACTORY_H_

#include "openmm/KernelFactory.h"

namespace OpenMM {

/**
 * This KernelFactory creates kernels for the reference implementation of the ATMMetaForce plugin.
 */

class ReferenceATMMetaForceKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEATMMETAFORCEKERNELFACTORY_H_*/
