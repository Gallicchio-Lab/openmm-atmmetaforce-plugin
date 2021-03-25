%module atmmetaforce

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "ATMMetaForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
from simtk.unit import *
import math
%}

%include "ATMMetaForceUtils.py"

/* 
 * Add units to function outputs.
*/
%pythonappend ATMMetaForcePlugin::ATMMetaForce::getParticleParameters(int index, int &particle, int& dx, int& dy, int& dz) const %{
    val[1] = Quantity(val[1], nanometer)
    val[2] = Quantity(val[2], nanometer)
    val[3] = Quantity(val[3], nanometer)
%}

%pythonappend ATMMetaForcePlugin::ATMMetaForce::getPerturbationEnergy(OpenMM::Context& context) const %{
   val=Quantity(val, kilojoules_per_mole)
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace ATMMetaForcePlugin {

class ATMMetaForce : public OpenMM::Force {
public:
    ATMMetaForce();

    int getNumParticles() const;
    int addParticle(int particle, double dx, double dy, double dz);
    void setParticleParameters(int index, int particle, double dx, double dy, double dz);
    void updateParametersInContext(OpenMM::Context& context);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply int& OUTPUT {int& particle};
    %apply double& OUTPUT {double& dx};
    %apply double& OUTPUT {double& dy};
    %apply double& OUTPUT {double& dz};
    void getParticleParameters(int index, int& particle, double& dx, double &dy, double &dz) const;
    %clear int& particle;
    %clear double& dx;
    %clear double& dy;
    %clear double& dz;

    double getPerturbationEnergy(OpenMM::Context& context) const;

    /**
     * get/set methods for the softplus alchemical potential function 
     * TODO: add units to alpha, U0, and W0
     */
    void setLambda1(double lambda1_t);
    double getLambda1() const;
    void setLambda2(double lambda2_t);
    double getLambda2() const;
    void setAlpha(double alpha_t);
    double getAlpha() const;
    void setU0(double u0_t);
    double getU0() const;
    void setW0(double w0_t);
    double getW0() const ;

    /**
     * get/set methods for the soft-core parameters
     * TODO: add units to Umax and Ubcore
     */
    double getUmax(void) const;
    void setUmax(double um);
    double getAcore(void) const;
    void setAcore(double a);
    double getUbcore(void) const;
    void setUbcore(double ub);    

    /*
     * Add methods for casting a Force to an ATMMetaForce.
    */
    %extend {
        static ATMMetaForcePlugin::ATMMetaForce& cast(OpenMM::Force& force) {
            return dynamic_cast<ATMMetaForcePlugin::ATMMetaForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<ATMMetaForcePlugin::ATMMetaForce*>(&force) != NULL);
        }
    }
};

}
