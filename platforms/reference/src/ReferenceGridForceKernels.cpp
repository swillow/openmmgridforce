/* -------------------------------------------------------------------------- *
 *                               OpenMMGridForce                              *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "ReferenceGridForceKernels.h"
#include "GridForce.h"

#include "openmm/OpenMMException.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/ReferencePlatform.h"

#include <cmath>

using namespace OpenMM;
using namespace std;

namespace GridForcePlugin {

// The length unit is nm
static vector<Vec3> &extractPositions(ContextImpl &context) {
    ReferencePlatform::PlatformData *data = reinterpret_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
    return *((vector<Vec3> *)data->positions);
}

static vector<Vec3> &extractForces(ContextImpl &context) {
    ReferencePlatform::PlatformData *data = reinterpret_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
    return *((vector<Vec3> *)data->forces);
}

/*
    OpenMM Grid Force
*/
void ReferenceCalcGridForceKernel::initialize(const System &system,
                                              const GridForce &grid_force) {
    // Initialize Nonbond parameters.
    grid_force.getGridParameters(g_counts, g_spacing, g_vals, g_scaling_factors);
}

double ReferenceCalcGridForceKernel::execute(ContextImpl &context,
                                             bool includeForces,
                                             bool includeEnergy) {
    
    vector<Vec3> &posData = extractPositions(context);
    vector<Vec3> &forceData = extractForces(context);
    
    const int nyz = g_counts[1] * g_counts[2];
    Vec3 hCorner(g_spacing[0] * (g_counts[0] - 1),
                 g_spacing[1] * (g_counts[1] - 1),
                 g_spacing[2] * (g_counts[2] - 1));

    double energy = 0.0;
    int natom_lig = g_scaling_factors.size();

    for (int ia = 0; ia < natom_lig; ++ia) {
        if (g_scaling_factors[ia] == 0.0) continue;

        Vec3 pi = posData[ia];
        bool is_inside = true;
        for (int k = 0; k < 3; ++k) {
            if (pi[k] > 0.0 && pi[k] < hCorner[k])
                continue;
            else
                is_inside = false;
        }

        if (is_inside) {
            int ix = (int)(pi[0] / g_spacing[0]);
            int iy = (int)(pi[1] / g_spacing[1]);
            int iz = (int)(pi[2] / g_spacing[2]);

            int im = ix * nyz + iy * g_counts[2] + iz;
            int imp = im + g_counts[2];  // iy --> iy + 1
            int ip = im + nyz;           // (ix --> ix+1)
            int ipp = ip + g_counts[2];  // (ix, iy) --> (ix+1, iy+1)

            // Corners of the box surrounding the point
            double vmmm = g_vals[im];
            double vmmp = g_vals[im + 1];   // iz --> iz+1
            double vmpm = g_vals[imp];      // iy --> iy + 1
            double vmpp = g_vals[imp + 1];  // (iy,iz)-->(iy+1, iz+1)

            double vpmm = g_vals[ip];
            double vpmp = g_vals[ip + 1];
            double vppm = g_vals[ipp];
            double vppp = g_vals[ipp + 1];

            // Fraction within the box
            double fx = (pi[0] - ix * g_spacing[0]) / g_spacing[0];
            double fy = (pi[1] - iy * g_spacing[1]) / g_spacing[1];
            double fz = (pi[2] - iz * g_spacing[2]) / g_spacing[2];

            // Fraction ahead
            double ax = 1.0 - fx;
            double ay = 1.0 - fy;
            double az = 1.0 - fz;

            // Trillinear interpolation for energy
            double vmm = az * vmmm + fz * vmmp;
            double vmp = az * vmpm + fz * vmpp;
            double vpm = az * vpmm + fz * vpmp;
            double vpp = az * vppm + fz * vppp;

            double vm = ay * vmm + fy * vmp;
            double vp = ay * vpm + fy * vpp;

	        double enr = g_scaling_factors[ia] * (ax * vm + fx * vp);

            energy += enr;

            // x coordinate
            double dvdx = -vm + vp;
            // y coordinate
            double dvdy = (-vmm + vmp) * ax + (-vpm + vpp) * fx;
            // z coordinate
            double dvdz = ((-vmmm + vmmp) * ay + (-vmpm + vmpp) * fy) * ax +
                          ((-vpmm + vpmp) * ay + (-vppm + vppp) * fy) * fx;
            Vec3 grd(dvdx / g_spacing[0], dvdy / g_spacing[1], dvdz / g_spacing[2]);
            forceData[ia] -= g_scaling_factors[ia] * grd;
        } else {
            double kval = 10000.0;  // kJ/mol nm**2
            Vec3 grd(0.0, 0.0, 0.0);
            for (int k = 0; k < 3; k++) {
                double dev = 0.0;
                if (pi[k] < 0.0) {
                    dev = pi[k];
                } else if (pi[k] > hCorner[k]) {
                    dev = pi[k] - hCorner[k];
                }
                energy += 0.5 * kval * dev * dev;
                grd[k] = kval * dev;
            }

            forceData[ia] -= g_scaling_factors[ia] * grd;
        }
    }

    return static_cast<double>(energy);
}

void ReferenceCalcGridForceKernel::copyParametersToContext(ContextImpl &context,
                                                           const GridForce &grid_force) {
    grid_force.getGridParameters(g_counts, g_spacing, g_vals, g_scaling_factors);
}



}  // namespace AlGDockPlugin
