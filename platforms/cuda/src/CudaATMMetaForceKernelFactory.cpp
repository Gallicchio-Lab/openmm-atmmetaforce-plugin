#include <exception>

#include "CudaATMMetaForceKernelFactory.h"
#include "CommonATMMetaForceKernels.h"
#include "CudaATMMetaForceKernels.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace ATMMetaForcePlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        CudaATMMetaForceKernelFactory* factory = new CudaATMMetaForceKernelFactory();
        platform.registerKernelFactory(CalcATMMetaForceKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerATMMetaForceCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* CudaATMMetaForceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
  CudaContext& cu = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
  if (name == CalcATMMetaForceKernel::Name())
    return new CudaCalcATMMetaForceKernel(name, platform, cu);
  throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
