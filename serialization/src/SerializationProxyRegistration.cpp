#include "openmm/serialization/SerializationProxy.h"

#include "ATMMetaForceProxy.h"
#include "ATMMetaForce.h"

#if defined(WIN32)
#include <windows.h>
    extern "C" void registerSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerSerializationProxies();
        return TRUE;
    }
#else
extern "C" void __attribute__((constructor)) registerSerializationProxies();
#endif

using namespace ATMMetaForcePlugin;

extern "C" void registerSerializationProxies() {
    OpenMM::SerializationProxy::registerProxy(typeid(ATMMetaForce), new ATMMetaForceProxy());
}
