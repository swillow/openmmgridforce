#ifndef OPENMM_GRIDFORCE_IMPL_H_
#define OPENMM_GRIDFORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                                OpenMMGridForce                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
 * Contributors:                                                              *
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



#include "GridForce.h"
#include "openmm/Kernel.h"
#include "openmm/internal/ForceImpl.h"
#include <string>
#include <utility>
#include <vector>


namespace GridForcePlugin {

/**
 * This is the internal implementation of AlGDockGridForce.
 */

class OPENMM_EXPORT_GRIDFORCE GridForceImpl : public OpenMM::ForceImpl {
   public:
    GridForceImpl(const GridForce &owner);
    ~GridForceImpl();
    void initialize(OpenMM::ContextImpl &context);
    const GridForce &getOwner() const {
        return owner;
    }
    void updateContextState(OpenMM::ContextImpl &context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(OpenMM::ContextImpl &context, 
                                bool includeForces, 
                                bool includeEnergy, 
                                int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>();  // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();

    void updateParametersInContext(OpenMM::ContextImpl &context);

   private:
    const GridForce &owner;
    OpenMM::Kernel kernel;
};

}  // namespace GridForcePlugin

#endif /*OPENMM_GRIDFORCE_IMPL_H_*/
