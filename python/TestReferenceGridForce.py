from __future__ import print_function

import unittest
from openmm import app
import openmm as omm
from openmm import unit
from openmm.vec3 import Vec3
import sys
import gridforceplugin
import numpy as np
#import algdock


def grid_read(FN):
    """
    Reads a grid in a netcdf format
    The multiplier affects the origin and spacing
    """

    if FN is None:
        raise Exception('File is not defined')
    elif FN.endswith('.nc'):
        from netCDF4 import Dataset
        grid_nc = Dataset(FN, 'r')
        data = {}
        for key in list(grid_nc.variables):
            data[key] = np.array(grid_nc.variables[key][:][0][:])
        grid_nc.close()
    else:
        raise Exception('File type not supported')

    return data


def getGridForce(FN, unit_conversion):

    data = grid_read(FN)

    force = gridforceplugin.GridForce()
    nx = int (data['counts'][0])
    ny = int (data['counts'][1])
    nz = int (data['counts'][2])
    force.addGridCounts(nx, ny, nz)
    # data['spacing] length unit is A
    # A -> nm: 0.1
    data['spacing'] *= 0.1
    force.addGridSpacing(data['spacing'][0],
                     data['spacing'][1],
                     data['spacing'][2])
    # vals:
    data['vals'] *= unit_conversion
    for val in data['vals']:
        force.addGridValue(val)

    return force


class TestCustomForce (unittest.TestCase):
    """This tests the Reference implementation of ReferenceAlGDockSphereForce."""

    def testGrid(self):
        prmtop = app.AmberPrmtopFile('prmtopcrd/ligand.prmtop')
        inpcrd = app.AmberInpcrdFile('prmtopcrd/ligand.trans.inpcrd')
        system = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds,
                                     implicitSolvent=None)

        # val unit = kcal/mol/e --> kJ/mol/e
        unit_conversion = 4.184
        force = getGridForce('grids/direct_ele.nc', unit_conversion)

        for chg in prmtop._prmtop.getCharges():
            force.addScalingFactor(chg)
        system.addForce(force)

        # val unit = sqrt(kcal/mol)/A^6 --> sqrt(kJ/mol)/nm^6
        unit_conversion = np.sqrt(4.184)*1.0e6
        force = getGridForce('grids/LJr.nc', unit_conversion)

        for rVdw, eps in prmtop._prmtop.getNonbondTerms():
            # rVdw : nm
            # eps  : kJ/mol
            ljr_scale = np.sqrt(eps)*(2.0*rVdw)**6
            force.addScalingFactor(ljr_scale)
        system.addForce(force)

        # val unit = sqrt(kcal/mol)/A^3 --> sqrt(kJ/mol)/nm^3
        unit_conversion = np.sqrt(4.184)*1.0e3
        force = getGridForce('grids/LJa.nc', unit_conversion)

        for rVdw, eps in prmtop._prmtop.getNonbondTerms():
            # rVdw : nm
            # eps  : kJ/mol
            lja_scale = np.sqrt(eps)*(2.0*rVdw)**3
            force.addScalingFactor(lja_scale)
        system.addForce(force)

        integrator = omm.LangevinIntegrator(300.0*unit.kelvin,
                                            1.0/unit.picoseconds,
                                            2.0*unit.femtoseconds)
        platform = omm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(
            prmtop.topology, system, integrator, platform)
        simulation.context.setPositions(inpcrd.positions)
        state = simulation.context.getState(
            getForces=True, getEnergy=True, getPositions=False)
        potential_energy = state.getPotentialEnergy()

        print('potential_energy ', potential_energy)


if __name__ == '__main__':
    unittest.main()
