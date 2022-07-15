#ifndef ATMMETAFORCE_HIPKERNELS_H_
#define ATMMETAFORCE_HIPKERNELS_H_

#include "ATMMetaForceKernels.h"
#include "CommonATMMetaForceKernels.h"
#include "openmm/hip/HipPlatform.h"
#include "openmm/hip/HipContext.h"

using namespace OpenMM;

namespace ATMMetaForcePlugin {
  
class HipCalcATMMetaForceKernel : public CommonCalcATMMetaForceKernel {
public:
    HipCalcATMMetaForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMMetaForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<HipPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

}

#endif /* ATMMETAFORCE_HIPKERNELS_H_ */
