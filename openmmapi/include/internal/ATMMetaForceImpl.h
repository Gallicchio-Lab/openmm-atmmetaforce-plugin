#ifndef OPENMM_ATMMETAFORCEFORCEIMPL_H_
#define OPENMM_ATMMETAFORCEFORCEIMPL_H_


#include "ATMMetaForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "internal/windowsExportATMMetaForce.h"
#include <utility>
#include <set>
#include <string>

namespace ATMMetaForcePlugin {

class System;

/**
 * This is the internal implementation of ATMMetaForce.
 */

//class OPENMM_EXPORT_ATMMETAFORCE ATMMetaForceImpl : public OpenMM::ForceImpl {
class OPENMM_EXPORT_ATMMETAFORCE ATMMetaForceImpl : public OpenMM::ForceImpl {
public:
    ATMMetaForceImpl(const ATMMetaForce& owner);
    ~ATMMetaForceImpl();
    void initialize(OpenMM::ContextImpl& context);
    const ATMMetaForce& getOwner() const {
        return owner;
    }
    void updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(OpenMM::ContextImpl& context);
    double getPerturbationEnergy() const {
       return PerturbationEnergy;
    }
private:
    const ATMMetaForce& owner;
    OpenMM::Kernel kernel;
    OpenMM::System innerSystem1, innerSystem2;
    OpenMM::VerletIntegrator innerIntegrator1, innerIntegrator2;
    OpenMM::Context *innerContext1, *innerContext2;
    double PerturbationEnergy;
};

} // namespace ATMMetaPlugin

#endif /*OPENMM_ATMMETAFORCEIMPL_H_*/
