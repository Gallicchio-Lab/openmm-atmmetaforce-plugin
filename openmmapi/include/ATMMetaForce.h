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
#include "internal/windowsExportATMMetaForce.h"

namespace ATMMetaForcePlugin {

/**
 * This class implements ...
 */

class OPENMM_EXPORT_ATMMETAFORCE ATMMetaForce : public OpenMM::Force {
public:
    /**
     * Create an ATMMetaForce.
     */
    ATMMetaForce();
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
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
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
      return false; //TODO
    }

  /**
   * returns the perturbation energy
   */
    double getPerturbationEnergy(const OpenMM::Context& context) const;


    /**
     * get/set methods for the softplus alchemical potential function parameters 
     */
    void setLambda1(double lambda1_t) {
      lambda1 = lambda1_t;
    }
    double getLambda1() const {
      return lambda1;
    }
    void setLambda2(double lambda2_t) {
      lambda2 = lambda2_t;
    }
    double getLambda2() const {
      return lambda2;
    }
    void setAlpha(double alpha_t) { //in 1/kjmol
      alpha = alpha_t;
    }
    double getAlpha() const {
      return alpha;
    }
    void setU0(double u0_t) { //in kjmol
      u0 = u0_t;
    }
    double getU0() const {
      return u0;
    }
    void setW0(double w0_t){
      w0 = w0_t;
    }
    double getW0() const {
      return w0;
    }

    /**
     * get/set methods for the soft-core parameters
     */
    double getUmax(void) const {
      return umax;
    }
    void setUmax(double um)  {
      umax = um;
    }
    double getAcore(void) const {
      return acore;
    }
    void setAcore(double a)  {
      acore = a;
    }
    double getUbcore(void) const {
      return ubcore;
    }
    void setUbcore(double ub)  {
      ubcore = ub;
    }
    
 protected:
  OpenMM::ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    std::vector<ParticleInfo> particles;

    //softplus parameters
    double lambda1, lambda2, alpha, u0, w0;
    //soft core parameters
    double umax, acore, ubcore;
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
