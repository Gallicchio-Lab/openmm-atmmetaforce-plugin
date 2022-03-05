#ifndef ATMMETAFORCE_CUDAKERNELS_H_
#define ATMMETAFORCE_CUDAKERNELS_H_

#include "ATMMetaForceKernels.h"
#include "CommonATMMetaForceKernels.h"
#include "openmm/cuda/CudaPlatform.h"
#include "openmm/cuda/CudaContext.h"

using namespace OpenMM;

namespace ATMMetaForcePlugin {
  
class CudaCalcATMMetaForceKernel : public CommonCalcATMMetaForceKernel {
public:
    CudaCalcATMMetaForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMMetaForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<CudaPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

}

#endif /* ATMMETAFORCE_CUDAKERNELS_H_ */
