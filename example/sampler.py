import simtk.unit as unit
from openmm import app
import openmm as omm
from openmm.vec3 import Vec3
from openmm.app import element as elem

import gridforceplugin

import numpy as np
import sys
import bat_mda as bat
import json
import os

rng_pcg64 = np.random.default_rng()


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
    nx = int(data['counts'][0])
    ny = int(data['counts'][1])
    nz = int(data['counts'][2])
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


def write_xyz(fxyz, istate, crd, en):

    title = "istate %10d %12.4f\n" % (istate, en)
    fxyz.write(title)
    natom = len(crd)
    line_natom = '%5d' % natom + '\n'
    fxyz.write(line_natom)
    for x, y, z in crd:
        line = '%12.7f' % x + '%12.7f' % y + '%12.7f' % z + '\n'
        fxyz.write(line)


class Sampler (object):

    def __init__(self,
                 system,
                 topology,
                 positions,
                 time_step):
        '''
        time_step : real (in femtosecond)
        '''
        self.do_replica_exchange = False
        self.do_genetic_MC = False

        self._system = system
        self._topology = topology
        self._positions = positions

        # --- Replica Exchange ---
        self._repX_temperatures = None
        self._repX_betas = None
        self._repX_interval = 0
        self._repX_simulations = None

        # --- Genetic MC
        self._zcoord = None
        self._primary_torsion = None

        integrator = omm.LangevinIntegrator(300.0*unit.kelvin,
                                            1.0/unit.picoseconds,
                                            time_step*unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)
        platform = omm.Platform.getPlatformByName('Reference')

        self._simulation = app.Simulation(
            topology, system, integrator, platform)
        self._simulation.context.setPositions(positions)
        

    def replica_exchange_init(self,
                              repX_temperatures=None,
                              repX_Interval=100,
                              reportInterval=1000,
                              reportFile=None):

        self.do_replica_exchange = True

        if repX_temperatures == None:
            print('Error in replica_exchange_init')
            sys.exit()

        self._repX_temperatures = repX_temperatures
        self._repX_betas = [1.0/(unit.MOLAR_GAS_CONSTANT_R*t)
                            for t in self._repX_temperatures]
        self._repX_interval = repX_Interval

        #
        self._repX_simulations = []
        for t in self._repX_temperatures:
            integrator = omm.LangevinIntegrator(t,
                                                1.0/unit.picoseconds,
                                                2.0*unit.femtoseconds)
            integrator.setConstraintTolerance(0.00001)
            platform = omm.Platform.getPlatformByName('Reference')
            simulation = app.Simulation(
                self._topology, self._system, integrator, platform)

            simulation.context.setPositions(self._positions)
            if reportFile != None:
                fname = reportFile + str(t._value) + '.dat'
                simulation.reporters.append(
                    app.StateDataReporter(fname,
                                          reportInterval,
                                          step=True,
                                          potentialEnergy=True,  # kJ/mole
                                          temperature=True,
                                          separator='     '))

            self._repX_simulations.append(simulation)

    def MD_with_step(self, nstep=0):

        if nstep > 0:
            for i, sim in enumerate(self._repX_simulations):
                sim.context.setVelocitiesToTemperature(
                    self._repX_temperatures[i])
                sim.step(nstep)  # MD job
        else:
            for i, sim in enumerate(self._repX_simulations):
                sim.context.setVelocitiesToTemperature(
                    self._repX_temperatures[i])
                sim.step(self._repX_interval)  # MD job

    def write_potential_energies(self, fout):

        for sim_i in self._repX_simulations:
            state_i = sim_i.context.getState(getEnergy=True)
            pot_i = state_i.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            fout.write("%12.4f" % pot_i)
        fout.write("\n")

    def write_xyz(self, fout):

        ###
        iend = len(self._repX_simulations)-1
        for istate in [0, iend]:
            sim = self._repX_simulations[istate]
            state = sim.context.getState(getPositions=True, getEnergy=True)
            pot = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            pos = state.getPositions(
                asNumpy=True).value_in_unit(unit.angstrom)
            write_xyz(fout, istate, pos, pot)

    def replica_exchange_temperature(self):

        if not self.do_replica_exchange:
            return

        numTemp = len(self._repX_temperatures)
        isel, jsel = rng_pcg64.integers(numTemp, size=2)

        if isel == jsel:
            jsel = isel+1 if (isel+1 < numTemp) else isel - 1

        state_i = self._repX_simulations[isel].context.getState(getEnergy=True,
                                                                getPositions=True)
        state_j = self._repX_simulations[jsel].context.getState(getEnergy=True,
                                                                getPositions=True)

        pot_i = state_i.getPotentialEnergy()
        pot_j = state_j.getPotentialEnergy()

# --- selection rule
#        log_ratio = -self._repX_betas[isel]*(pot_j-pot_i) \
#                    -self._repX_betas[jsel]*(pot_i-pot_j)
#
        log_ratio = (self._repX_betas[isel] -
                     self._repX_betas[jsel]) * (pot_i - pot_j)

        ratio = rng_pcg64.random()
        nacc = 0
        if log_ratio >= 0:
            nacc = 1
        elif ratio < np.exp(log_ratio):
            nacc = 1

        if nacc == 1:
            # the exchange is accepted
            xyz_i = state_i.getPositions()
            xyz_j = state_j.getPositions()
            self._repX_simulations[isel].context.setPositions(xyz_j)
            self._repX_simulations[jsel].context.setPositions(xyz_i)

        return nacc

    def genetic_MC_init(self):

        self.do_genetic_MC = True

        natom = len(self._positions)
        ibonds = np.array([[b[0].index, b[1].index]
                           for b in self._topology.bonds()])
        masses = np.zeros(natom, dtype=np.float)

        for ia in range(natom):
            masses[ia] = self._system.getParticleMass(ia)._value

        self._zcoord, self._primary_torsion = \
            bat.get_zcoord(masses, ibonds)

    def genetic_MC_mutation(self):
        if not self.do_genetic_MC:
            return

        numConf = len(self._repX_temperatures)  # number of Conf.
        isel, jsel = rng_pcg64.integers(numConf, size=2)

        if isel == jsel:
            jsel = isel+1 if (isel+1 < numConf) else isel-1

        # isel is a system with lower temperature
        if jsel < isel:
            (isel, jsel) = (jsel, isel)

        state_i = self._repX_simulations[isel].context.getState(getEnergy=True,
                                                                getPositions=True)
        state_j = self._repX_simulations[jsel].context.getState(getEnergy=True,
                                                                getPositions=True)
        pot_i = state_i.getPotentialEnergy()
        pot_j = state_j.getPotentialEnergy()

        # --- Get BAT coordinats from XYZ
        xyz_i = state_i.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        xyz_j = state_j.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        # bat_i : np.array
        bat_i = bat.get_bat_from_xyz(
            xyz_i, self._zcoord, self._primary_torsion)
        bat_j = bat.get_bat_from_xyz(
            xyz_j, self._zcoord, self._primary_torsion)

        # --- Do Crossover using BAT coordinates
        nzcrd = len(self._zcoord)
        tor_i = bat_i[2*nzcrd+9:]
        tor_j = bat_j[2*nzcrd+9:]

        tor_a = tor_i
        tor_b = tor_j

        """
        unique_primary_torsion_indices = list(set(self._primary_torsion))
        icut = rng_pcg64.integers(
            len(unique_primary_torsion_indices), size=1)[0]
        icut2 = unique_primary_torsion_indices[icut]

        tor_a[icut2] = tor_j[icut2]
        tor_b[icut2] = tor_i[icut2]
        """
        icut = rng_pcg64.integers(nzcrd, size=1)[0]
        tor_a[icut] = tor_j[icut]
        bat_i[2*nzcrd+9:] = tor_a

        # --- Get new XYZ coordinates from new BATs
        xyz_a = unit.Quantity(bat.get_xyz_from_bat(bat_i, self._zcoord, self._primary_torsion),
                              unit.angstrom)


        # --- Estimate new potential energeis using new XYZs
        self._simulation.context.setPositions(xyz_a)
        state_a = self._simulation.context.getState(getEnergy=True)
        pot_a = state_a.getPotentialEnergy()

        # --- selection rule
        log_ratio = -self._repX_betas[isel]*(pot_a-pot_i)


        ratio = rng_pcg64.random()
        nacc = 0
        if log_ratio >= 0:
            nacc = 1
        elif ratio < np.exp(log_ratio):
            nacc = 1

        if log_ratio >= 50:
            print("I THINK ERROR mutation", pot_a,
                  pot_i, icut, log_ratio, isel, jsel)
            nacc = 0

        if nacc == 1:
            # accepted
            #print ('isel', isel, jsel, 'icut ', icut, ndof+offset, 'E_before', pot_i, pot_j, 'E_after', pot_a, pot_b)
            self._repX_simulations[isel].context.setPositions(xyz_a)
 

        return nacc

    def genetic_MC_crossover(self):

        if not self.do_genetic_MC:
            return

        numConf = len(self._repX_temperatures)  # number of Conf.
        isel, jsel = rng_pcg64.integers(numConf, size=2)

        if isel == jsel:
            jsel = isel+1 if (isel+1 < numConf) else isel-1

        # isel is a system with lower temperature
        if jsel < isel:
            (isel, jsel) = (jsel, isel)

        state_i = self._repX_simulations[isel].context.getState(getEnergy=True,
                                                                getPositions=True)
        state_j = self._repX_simulations[jsel].context.getState(getEnergy=True,
                                                                getPositions=True)
        pot_i = state_i.getPotentialEnergy()
        pot_j = state_j.getPotentialEnergy()

        # --- Get BAT coordinats from XYZ
        xyz_i = state_i.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        xyz_j = state_j.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        # bat_i : np.array
        bat_i = bat.get_bat_from_xyz(
            xyz_i, self._zcoord, self._primary_torsion)
        bat_j = bat.get_bat_from_xyz(
            xyz_j, self._zcoord, self._primary_torsion)

        # --- Do Crossover using BAT coordinates
        nzcrd = len(self._zcoord)
        tor_i = bat_i[2*nzcrd+9:]
        tor_j = bat_j[2*nzcrd+9:]

        icut = rng_pcg64.integers(nzcrd, size=1)[0]
        """
        unique_primary_torsion_indices = list(set(self._primary_torsion))
        icut = rng_pcg64.integers(
            len(unique_primary_torsion_indices), size=1)[0]
        icut2 = unique_primary_torsion_indices[icut]

        tor_a = np.array(list(tor_i[:icut2]) + list(tor_j[icut2:]))
        tor_b = np.array(list(tor_j[:icut2]) + list(tor_i[icut2:]))
        """
        bat_i[2*nzcrd+9:] = np.array(list(tor_i[:icut])+list(tor_j[icut:]))
        #bat_j[2*nzcrd+9:] = tor_b

        # --- Get new XYZ coordinates from new BATs
        xyz_a = unit.Quantity(bat.get_xyz_from_bat(bat_i, self._zcoord, self._primary_torsion),
                              unit.angstrom)
        # xyz_b = unit.Quantity(bat.get_xyz_from_bat(bat_j, self._zcoord, self._primary_torsion),
        #                      unit.angstrom)

        # --- Estimate new potential energeis using new XYZs
        self._simulation.context.setPositions(xyz_a)
        state_a = self._simulation.context.getState(getEnergy=True)
        pot_a = state_a.getPotentialEnergy()

        # self._simulation.context.setPositions(xyz_b)
        #state_b = self._simulation.context.getState(getEnergy=True)
        #pot_b = state_b.getPotentialEnergy()

        # if pot_a > pot_b:
        #    (pot_a, pot_b) = (pot_b, pot_a)
        #    (xyz_a, xyz_b) = (xyz_b, xyz_a)

        # --- selection rule
        log_ratio = -self._repX_betas[isel]*(pot_a-pot_i)
#                    - self._repX_betas[jsel]*(pot_b-pot_j)

        nacc = 0
        if log_ratio >= 0:
            nacc = 1
        else:
            ratio = rng_pcg64.random()

            if ratio < np.exp(log_ratio):
                nacc = 1

        if log_ratio >= 30:
            print("I THINK ERROR cross", pot_a,
                  pot_i, icut, log_ratio, isel, jsel)
            nacc = 0

        if nacc == 1:
            # accepted
            #print ('isel', isel, jsel, 'icut ', icut, ndof+offset, 'E_before', pot_i, pot_j, 'E_after', pot_a, pot_b)
            self._repX_simulations[isel].context.setPositions(xyz_a)
#            self._repX_simulations[jsel].context.setPositions(xyz_b)
#        else:
#            self._repX_simulations[isel].context.setPositions(
#                unit.Quantity(xyz_i, unit.angstrom))
#           self._repX_simulations[jsel].context.setPositions(
#                unit.Quantity(xyz_j, unit.angstrom))

        return nacc



def main_run(inp_fname):



    # Read input.json
    with open(inp_fname) as f:
        data = json.load(f)

    run_job = data['run_job']
    nstate = data['nstate']
    ntrial_repX = data['ntrial_repX']
    ntrial_gMC = data['ntrial_gMC']

    nstep_equil = data['nstep_equil']

    work_dir = data['work_dir'] + '/' + run_job + '/' + \
        str(nstate) + '_' + str(ntrial_repX) + '_' + str(ntrial_gMC) + '/'
    if not os.path.isdir(work_dir):
        cmd = 'mkdir -p ' + work_dir
        os.system(cmd)

    solvent = None
    # if data[run_job]['solvation'] == 'Desolvated':
    #    solvent = app.OBC2

    # ('prmtopcrd/ligand.prmtop')
    prmtop = app.AmberPrmtopFile(data['dir']['ligand_prmtop'])
    # ('prmtopcrd/ligand.inpcrd')
    inpcrd = app.AmberInpcrdFile(data['dir']['ligand_inpcrd'])

    masses = prmtop._prmtop.getMasses()
    pos = inpcrd.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

    # Get Center of Mass of Ligand
    R0 = np.array([0.0, 0.0, 0.0])
    for ia in range(len(masses)):
        R0 += masses[ia]*pos[ia]

    R0 /= sum(masses)

    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds,
                                 implicitSolvent=solvent)

    # Heavy Hydrogen Mass
    Hmass = data[run_job]['H_mass']
    delta_t = data[run_job]['delta_t']

    if Hmass:
        for atom in prmtop.topology.atoms():
            if atom.element == elem.hydrogen:
                system.setParticleMass(atom.index, Hmass)
    
    if run_job == 'CD':
        # R0 = data['CD']['site_center']  # length unit is nm
        #Rcenter = Vec3(R0[0], R0[1], R0[2])
        #Rmax = data['CD']['site_max_R']
        #system.addForce(algdockplugin.AlGDockSphereForce(Rcenter, Rmax))

        # Grid: direct_elec
        unit_conversion = 4.184  # kcal/mol/e --> kJ/mol/e
        force = getGridForce(
            data['grids']['direct_elec'], unit_conversion)

        for chg in prmtop._prmtop.getCharges():
            force.addScalingFactor(chg)
        system.addForce(force)

        # sqrt(kcal/mol)/A^6 --> sqrt(kJ/mol)/nm^6
        unit_conversion = np.sqrt(4.184)*1.0e6
        force = getGridForce(
            data['grids']['LJr'], unit_conversion)

        for rVdw, eps in prmtop._prmtop.getNonbondTerms():
            #rVdw : nm
            # eps : kJ/mol
            ljr_scale = np.sqrt(eps)*(2.0*rVdw)**6
            force.addScalingFactor(ljr_scale)
        system.addForce(force)

        # sqrt(kcal/mol)/A^3 --> sqrt(kJ/mol)/nm^3
        unit_conversion = np.sqrt(4.184)*1.0e3
        force = getGridForce(
            data['grids']['LJa'], unit_conversion)

        for rVdw, eps in prmtop._prmtop.getNonbondTerms():
            #rVdw : nm
            # eps : kJ/mol
            lja_scale = np.sqrt(eps)*(2.0*rVdw)**3
            force.addScalingFactor(lja_scale)
        system.addForce(force)

    # traj = mdtraj.load('prmtopcrd/ligand.inpcrd',
    #                   top='prmtopcrd/ligand.prmtop')

    sampler = Sampler(system,
                      prmtop.topology,
                      inpcrd.positions,
                      delta_t)

    reportFile = None  # work_dir+'repX'
    logFile = work_dir+'sample.log'
    potFile = work_dir+'energy.log'
    xyzFile = work_dir+'xyz.inpcrd'

    minT = 300.0*unit.kelvin
    maxT = 600.0*unit.kelvin
    numT = nstate
    nstep_md = 100  # ReplicaExchangeInterval

    repX_temperatures = [minT + t*(maxT-minT)/(numT-1) for t in range(numT)]

    sampler.replica_exchange_init(repX_temperatures,
                                  repX_Interval=nstep_md,
                                  reportFile=reportFile)

    if ntrial_gMC > 0:
        sampler.genetic_MC_init()

    # Equilibrium
    sampler.MD_with_step(nstep_equil)

    # Production
    nacc = 0
    nacc_gmc = 0

    fout_pot = open(potFile, 'w', 1)
    fout = open(logFile, 'w', 1)
    fout_xyz = open(xyzFile, 'w', 1)

    ntrial = 10000
    nstep_md = data['nstep_MD'] - 2*ntrial_gMC
    for i in range(ntrial):

        # (1) Do replica exchange
        for it in range(ntrial_repX):
            nacc += sampler.replica_exchange_temperature()

        # (2) Do Genetic MC
        for it in range(ntrial_gMC):
            nacc_gmc += sampler.genetic_MC_crossover()
            nacc_gmc += sampler.genetic_MC_mutation()

        # (3) Perform MD simulations (MDstep) with N different thermal states
        sampler.MD_with_step(nstep_md)  # 100 MDstep * 4 fs

        sampler.write_potential_energies(fout_pot)
        if (i+1) % 10 == 0:
            sampler.write_xyz(fout_xyz)

        if (i+1) % 100 == 0:
            ratio_repX = 0.0
            ratio_gMC = 0.0
            if ntrial_repX > 0:
                ratio_repX = float(nacc)/((i+1)*ntrial_repX)
            if ntrial_gMC > 0:
                ratio_gMC = float(nacc_gmc) / ((i+1)*ntrial_gMC)
            fout.write('acceptance_at_ntrial: %8d %12.4f %12.4f \n' %
                       (i+1, ratio_repX, ratio_gMC))

    fout.close()
    fout_pot.close()
    fout_xyz.close()


#        if nacc_gmc > 10:
#            traj.save ('test.xyz')
#            sys.exit()


#    print ('acceptance ', float(nacc)/ntrial)
#    print ('acceptance ', float(nacc_gmc)/ntrial)
if __name__ == '__main__':
    import getopt
    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "i:", ["ifile="])

    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inp_fname = arg

    main_run(inp_fname)
