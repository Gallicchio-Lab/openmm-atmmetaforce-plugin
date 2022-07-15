
#include <exception>

#include "HipATMMetaForceKernelFactory.h"
#include "CommonATMMetaForceKernels.h"
#include "HipATMMetaForceKernels.h"
#include "openmm/hip/HipContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace ATMMetaForcePlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("HIP");
        HipATMMetaForceKernelFactory* factory = new HipATMMetaForceKernelFactory();
        platform.registerKernelFactory(CalcATMMetaForceKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerATMMetaForceHipKernelFactories() {
    try {
        Platform::getPlatformByName("HIP");
    }
    catch (...) {
        Platform::registerPlatform(new HipPlatform());
    }
    registerKernelFactories();
}

KernelImpl* HipATMMetaForceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    HipContext& cl = *static_cast<HipPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == CalcATMMetaForceKernel::Name())
        return new HipCalcATMMetaForceKernel(name, platform, cl);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
