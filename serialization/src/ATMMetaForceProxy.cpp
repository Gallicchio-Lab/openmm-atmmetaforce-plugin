#include "ATMMetaForceProxy.h"
#include "ATMMetaForce.h"

#include <vector>

using namespace ATMMetaForcePlugin;

ATMMetaForceProxy::ATMMetaForceProxy() : OpenMM::SerializationProxy("ATMMetaForce") {
}

void ATMMetaForceProxy::serialize(const void* object, OpenMM::SerializationNode& node) const {
    node.setIntProperty("version", 0);
    const ATMMetaForce& force = *reinterpret_cast<const ATMMetaForce*>(object);

    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("name", force.getName());

    node.setDoubleProperty("lambda1", force.getDefaultLambda1());
    node.setDoubleProperty("lambda2", force.getDefaultLambda2());
    node.setDoubleProperty("alpha", force.getDefaultAlpha());
    node.setDoubleProperty("u0", force.getDefaultU0());
    node.setDoubleProperty("w0", force.getDefaultW0());
    node.setDoubleProperty("uMax", force.getDefaultUmax());
    node.setDoubleProperty("ubCore", force.getDefaultUbcore());
    node.setDoubleProperty("aCore", force.getDefaultAcore());
    node.setDoubleProperty("direction", force.getDefaultDirection());

    OpenMM::SerializationNode& variableForceGroups = node.createChildNode("GlobalParameters");
    for (const auto i : force.getVariableForceGroups()) {
        variableForceGroups.createChildNode("Parameter").setIntProperty("group", i);
    }

    OpenMM::SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int particle;
        double dx, dy, dz;
        force.getParticleParameters(i, particle, dx, dy, dz);
        particles.createChildNode("Particle").setIntProperty("particle", particle).setDoubleProperty("dx", dx).setDoubleProperty("dy", dy).setDoubleProperty("dz", dz);
    }
}

void* ATMMetaForceProxy::deserialize(const OpenMM::SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version != 0) throw OpenMM::OpenMMException("Unsupported version");

    ATMMetaForce* force = nullptr;
    try {
        std::vector<int> variableForceGroups;

        const OpenMM::SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            variableForceGroups.push_back(parameter.getIntProperty("group"));

        force = new ATMMetaForce(
            node.getDoubleProperty("lambda1"),
            node.getDoubleProperty("lambda2"),
            node.getDoubleProperty("alpha"),
            node.getDoubleProperty("u0"),
            node.getDoubleProperty("w0"),
            node.getDoubleProperty("uMax"),
            node.getDoubleProperty("ubCore"),
            node.getDoubleProperty("aCore"),
            node.getDoubleProperty("direction"),
            variableForceGroups
        );
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));

        const OpenMM::SerializationNode& particles = node.getChildNode("Particles");
        for (auto& particle : particles.getChildren())
            force->addParticle(particle.getIntProperty("particle"), particle.getDoubleProperty("dx"), particle.getDoubleProperty("dy"), particle.getDoubleProperty("dz"));

        return force;
    }
    catch (...) {
        if (force != nullptr)
            delete force;
        throw;
    }
}