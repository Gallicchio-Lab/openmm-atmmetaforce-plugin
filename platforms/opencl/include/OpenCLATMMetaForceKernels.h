#ifndef ATMMETAFORCE_OPENCLKERNELS_H_
#define ATMMETAFORCE_OPENCLKERNELS_H_

#include "ATMMetaForceKernels.h"
#include "CommonATMMetaForceKernels.h"
#include "openmm/opencl/OpenCLPlatform.h"
#include "openmm/opencl/OpenCLContext.h"

using namespace OpenMM;

namespace ATMMetaForcePlugin {
  
class OpenCLCalcATMMetaForceKernel : public CommonCalcATMMetaForceKernel {
public:
    OpenCLCalcATMMetaForceKernel(std::string name, const Platform& platform, ComputeContext& cc) : CommonCalcATMMetaForceKernel(name, platform, cc) {
    }
    ComputeContext& getInnerComputeContext(ContextImpl& innerContext) {
        return *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    }
};

}

#endif /* ATMMETAFORCE_OPENCLKERNELS_H_ */
