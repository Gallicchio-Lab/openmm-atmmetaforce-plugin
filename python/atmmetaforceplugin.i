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
    ATMMetaForce(double Lambda);
    ATMMetaForce(double Lambda1, double Lambda2, double Alpha, double U0, double W0, double Umax, double Ubcore, double Acore);
    ATMMetaForce(double Lambda1, double Lambda2, double Alpha, double U0, double W0, double Umax, double Ubcore, double Acore, double direction);

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

    static const std::string& Lambda1();
    static const std::string& Lambda2();
    static const std::string& Alpha();
    static const std::string& U0();
    static const std::string& W0();
    static const std::string& Umax();
    static const std::string& Ubcore();
    static const std::string& Acore();
    static const std::string& Direction();

    double getDefaultLambda1() const;
    double getDefaultLambda2() const;
    double getDefaultAlpha() const;
    double getDefaultU0() const;
    double getDefaultW0() const;
    double getDefaultUmax() const;
    double getDefaultUbcore() const;
    double getDefaultAcore() const;
    double getDefaultDirection() const;

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
