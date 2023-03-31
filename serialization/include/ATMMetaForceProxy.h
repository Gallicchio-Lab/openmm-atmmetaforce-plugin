#ifndef OPENMM_ATMMETAFORCE_PROXY_H_
#define OPENMM_ATMMETAFORCE_PROXY_H_

#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"

#include "internal/windowsExportATMMetaForce.h"

namespace ATMMetaForcePlugin {

/**
 * This is a proxy for serializing ATMMetaForce objects.
 */

class OPENMM_EXPORT_ATMMETAFORCE ATMMetaForceProxy : public OpenMM::SerializationProxy {
    public:
        ATMMetaForceProxy();
        void serialize(const void* object, OpenMM::SerializationNode& node) const override;
        void* deserialize(const OpenMM::SerializationNode& node) const override;
    };

} // namespace ATMMetaPlugin

#endif //OPENMM_ATMMETAFORCE_PROXY_H_
