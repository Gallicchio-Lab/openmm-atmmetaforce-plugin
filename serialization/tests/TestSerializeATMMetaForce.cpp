#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include "ATMMetaForce.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace ATMMetaForcePlugin;
using namespace OpenMM;  // needed for ASSERTS to work

void testSerialization() {
    const std::string expectedName = "MyATMMetaForce";
    const auto expectedForceGroup = 30;
    const auto expectedLambda1 = 0.0;
    const auto expectedLambda2 = 0.1;
    const auto expectedAlpha = 0.25;
    const auto expectedU0 = 0.5;
    const auto expectedW0 = 0.6;
    const auto expectedUmax = 200.0;
    const auto expectedUbcore = 100.0;
    const auto expectedAcore = 0.07;
    const auto expectedDirection=0.0;
    const std::vector<int> expectedVariableForceGroups={1};

    ATMMetaForce force = {
        expectedLambda1,
        expectedLambda2,
        expectedAlpha,
        expectedU0,
        expectedW0,
        expectedUmax,
        expectedUbcore,
        expectedAcore,
        expectedDirection,
        expectedVariableForceGroups
    };
    force.setForceGroup(expectedForceGroup);
    force.setName(expectedName);
    force.addParticle(0, 0.1, 0.2, 0.3);
    force.addParticle(2, 0.4, 0.5, 0.6);

    // Serialize and then deserialize it.

    std::stringstream buffer;
    XmlSerializer::serialize<ATMMetaForce>(&force, "Force", buffer);
    auto* copy = XmlSerializer::deserialize<ATMMetaForce>(buffer);

    ASSERT_EQUAL(expectedForceGroup, copy->getForceGroup());
    ASSERT_EQUAL(expectedName, copy->getName());

    ASSERT_EQUAL(expectedLambda1, copy->getDefaultLambda1());
    ASSERT_EQUAL(expectedLambda2, copy->getDefaultLambda2());
    ASSERT_EQUAL(expectedAlpha, copy->getDefaultAlpha());
    ASSERT_EQUAL(expectedU0, copy->getDefaultU0());
    ASSERT_EQUAL(expectedW0, copy->getDefaultW0());
    ASSERT_EQUAL(expectedUmax, copy->getDefaultUmax());
    ASSERT_EQUAL(expectedUbcore, copy->getDefaultUbcore());
    ASSERT_EQUAL(expectedAcore, copy->getDefaultAcore());
    ASSERT_EQUAL(expectedDirection, copy->getDefaultDirection());
    ASSERT(expectedVariableForceGroups == copy->getVariableForceGroups());

    ASSERT_EQUAL(copy->getNumParticles(), 2);
    int particle;
    double dx, dy, dz;;
    copy->getParticleParameters(0, particle, dx, dy, dz);
    ASSERT_EQUAL(0, particle);
    ASSERT_EQUAL(0.1, dx);
    ASSERT_EQUAL(0.2, dy);
    ASSERT_EQUAL(0.3, dz);
    copy->getParticleParameters(1, particle, dx, dy, dz);
    ASSERT_EQUAL(2, particle);
    ASSERT_EQUAL(0.4, dx);
    ASSERT_EQUAL(0.5, dy);
    ASSERT_EQUAL(0.6, dz);
}

int main() {
    try {
        testSerialization();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}