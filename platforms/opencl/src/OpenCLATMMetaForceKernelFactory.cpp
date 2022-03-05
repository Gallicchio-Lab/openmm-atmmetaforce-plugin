
#include <exception>

#include "OpenCLATMMetaForceKernelFactory.h"
#include "CommonATMMetaForceKernels.h"
#include "OpenCLATMMetaForceKernels.h"
#include "openmm/opencl/OpenCLContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace ATMMetaForcePlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("OpenCL");
        OpenCLATMMetaForceKernelFactory* factory = new OpenCLATMMetaForceKernelFactory();
        platform.registerKernelFactory(CalcATMMetaForceKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
    }
}

extern "C" OPENMM_EXPORT void registerATMMetaForceOpenCLKernelFactories() {
    try {
        Platform::getPlatformByName("OpenCL");
    }
    catch (...) {
        Platform::registerPlatform(new OpenCLPlatform());
    }
    registerKernelFactories();
}

KernelImpl* OpenCLATMMetaForceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    OpenCLContext& cl = *static_cast<OpenCLPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == CalcATMMetaForceKernel::Name())
        return new OpenCLCalcATMMetaForceKernel(name, platform, cl);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
