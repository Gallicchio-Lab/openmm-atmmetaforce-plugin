#ifndef OPENMM_ATMMETAFORCE_H_
#define OPENMM_ATMMETAFORCE_H_

/* -------------------------------------------------------------------------- *
 *                     OpenMM ATM Meta-Force Plugin                           *
 * -------------------------------------------------------------------------- *
 * This is a plugin of the OpenMM molecular simulation toolkit                *
 * (https://openmm.org) that implements the Alchemical Transfer Potential     *
 * for absolute and relative binding free energy estimation                   *
 * (https://arxiv.org/abs/2101.07894). The work is supported by a grant from  *
 * the National Science Foundation (NSF 1750511).                             *
 *                                                                            *
 * Portions copyright (c) 2021 by the Authors.                                *
 * Authors: Emilio Gallicchio                                                 *
 * Contributors: Peter Eastman                                                *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include <vector>
#include <string>
#include "internal/windowsExportATMMetaForce.h"
#include "ATMMetaForceVersion.h"

namespace ATMMetaForcePlugin {

/**
 * This class implements the Alchemical Transfer force
 */

class OPENMM_EXPORT_ATMMETAFORCE ATMMetaForce : public OpenMM::Force {
public:
    /**
     * Create an ATMMetaForce.
     */
 ATMMetaForce(double Lambda1, double Lambda2, double Alpha, double U0, double W0, double Umax, double Ubcore, double Acore, double direction, const std::vector<int>& VariableForceGroups) :
     defaultLambda1(Lambda1), defaultLambda2(Lambda2), defaultAlpha(Alpha), defaultU0(U0), defaultW0(W0),
       defaultUmax(Umax), defaultUbcore(Ubcore), defaultAcore(Acore), defaultDirection(direction), VariableForceGroups(VariableForceGroups) {
     }
    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add a particle to the force.
     *
     * @param particle    the index of the particle
     * @param dx, dy, dz  the displacement vector in nm
     * @return the index of the particle that was added
     */
    int addParticle(int particle, double dx, double dy, double dz);
    /**
     * Get the parameters for a particle
     * 
     * @param index      the index in the force for the particle for which to get parameters
     * @param particle   the index of the particle
     * @param dx, dy, dz the coordinates of the displacement vector in nm
     */
    void getParticleParameters(int index, int& particle, double& dx, double &dy, double &dz) const;
    /**
     * Set the force field parameters for a particle
     * 
     * @param index      the index in the force of the particle for which to set parameters
     * @param particle   the particle associated with this index
     * @param dx, dy, dz the coordinates of the displacement vector in nm
     */
    void setParticleParameters(int index, int particle, double dx, double dy, double dz);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method 
     * should be called after updating parameters with setParticleParameters() to copy them over to the Context.
     * The only information this method updates is the values of per-particle parameters.  The number of particles
     * cannot be changed.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
      return false; //the non-bonded force with PBC is in the system so it would be queried correctly
    }

  /**
   * returns the perturbation energy
   */
    double getPerturbationEnergy(const OpenMM::Context& context) const;


    /**
     * Names of the parameters
     */
    static const std::string& Lambda1() {
        static const std::string key = "ATMLambda1";
        return key;
    }
    static const std::string& Lambda2() {
        static const std::string key = "ATMLambda2";
        return key;
    }
    static const std::string& Alpha() {
        static const std::string key = "ATMAlpha";
        return key;
    }
    static const std::string& U0() {
        static const std::string key = "ATMU0";
        return key;
    }
    static const std::string& W0() {
        static const std::string key = "ATMW0";
        return key;
    }
    static const std::string& Umax() {
        static const std::string key = "ATMUmax";
        return key;
    }
    static const std::string& Ubcore() {
        static const std::string key = "ATMUbcore";
        return key;
    }
    static const std::string& Acore() {
        static const std::string key = "ATMAcore";
        return key;
    }
    static const std::string& Direction() {
        static const std::string key = "ATMDirection";
        return key;
    }

    static const std::string& Version() {
      //version id
      static const std::string version = ATMMETAFORCE_VERSION ;
      return version;
    }
    
    /**
     *  default values of the parameters
     */
    double getDefaultLambda1() const {
        return defaultLambda1;
    }
    double getDefaultLambda2() const {
        return defaultLambda2;
    }
    double getDefaultAlpha() const {
        return defaultAlpha;
    }
    double getDefaultU0() const {
        return defaultU0;
    }
    double getDefaultW0() const {
        return defaultW0;
    }
    double getDefaultUmax() const {
        return defaultUmax;
    }
    double getDefaultUbcore() const {
        return defaultUbcore;
    }
    double getDefaultAcore() const {
        return defaultAcore;
    }
    double getDefaultDirection() const {
        return defaultDirection;
    }

    const std::vector<int>& getVariableForceGroups() const {
      return VariableForceGroups;
    }

protected:
  OpenMM::ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    std::vector<ParticleInfo> particles;

    //softplus parameters
    double defaultLambda1, defaultLambda2, defaultAlpha, defaultU0, defaultW0;
    //soft core parameters
    double defaultUmax, defaultUbcore, defaultAcore;
    //alchemical direction parameter
    double defaultDirection;

    //the forces that are recalculated after the coordinate transformation
    std::vector<int> VariableForceGroups;

};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ATMMetaForce::ParticleInfo {
 public:
  int particle;
  double dx, dy, dz;
  ParticleInfo() {
    particle = -1;
    dx = dy = dz = 0.0;
  }
  ParticleInfo(int particle) : particle(particle) {
    dx = dy = dz = 0.0;
  }
  ParticleInfo(int particle, double dx, double dy, double dz) : particle(particle), dx(dx), dy(dy), dz(dz) {
  }
};
 
} // namespace ATMMetaForce

#endif /*OPENMM_ATMMETAFORCE_H_*/
