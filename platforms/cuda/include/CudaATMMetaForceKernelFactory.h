#ifndef OPENMM_CUDAATMMETAFORCEKERNELFACTORY_H_
#define OPENMM_CUDAATMMETAFORCEKERNELFACTORY_H_

#include "openmm/KernelFactory.h"

namespace OpenMM {

/**
 * This KernelFactory creates kernels for the CUDA implementation of the ATM Meta Force plugin.
 */

class CudaATMMetaForceKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAATMMETAFORCEKERNELFACTORY_H_*/
