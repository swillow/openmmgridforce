import numpy as np
import sys


def _distance(p1, p2):
    v1 = p2 - p1
    return np.sqrt(np.einsum('i,i', v1, v1))


def _angle(p1, p2, p3):
    v1 = p2 - p1
    v2 = p2 - p3
    val = np.einsum('i,i', v1, v2) / \
        np.sqrt(np.einsum('i,i', v1, v1)*np.einsum('i,i', v2, v2))
    return np.arccos(max(-1., min(1., val)))


def _dihedral(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    c1 = np.cross(b2, b3)
    c2 = np.cross(b1, b2)

    p1 = (b1*c1).sum(-1)
    p1 *= (b2*b2).sum(-1)**0.5
    p2 = (c1*c2).sum(-1)

    return np.arctan2(p1, p2)  # return in radian


def _sort_atoms_by_mass(atoms, reverse=False):

    return sorted(atoms, key=lambda a: (a[1], a[0]), reverse=reverse)


def sortFuncByMass(atom):
    return atom[1]


def _find_torsion(root, masses, bondToAtoms):

    torsions = []
    selected_atoms = list(root)
    while len(selected_atoms) < len(bondToAtoms):
        torsionAdded = False

        for a1 in selected_atoms:
            a0_list = _sort_atoms_by_mass([(a0, masses[a0])
                                           for a0 in bondToAtoms[a1]
                                           if (a0 not in selected_atoms)]
                                          )

            for a0, a0mass in a0_list:
                a2_list = _sort_atoms_by_mass([(a2, masses[a2]) for a2 in bondToAtoms[a1]
                                               if (a2 != a0) and
                                               (len(bondToAtoms[a2]) > 1) and
                                               (a2 in selected_atoms)])

                for a2, a2mass in a2_list:
                    a3_list = _sort_atoms_by_mass([(a3, masses[a3]) for a3 in bondToAtoms[a2]
                                                   if (a3 != a1) and
                                                   (a3 in selected_atoms)])

                    for a3, a3mass in a3_list:

                        torsions.append([a0, a1, a2, a3])
                        selected_atoms.append(a0)
                        torsionAdded = True
                        break  # out of the a3 loop
                    break  # out of the a2 loop

        if torsionAdded is False:
            print("Selected Atoms:", len(selected_atoms))
            selected_atoms.sort()
            print(selected_atoms)
            print("Torsions ", torsions)
            raise ValueError("Error")

    return np.array(torsions)


def get_zcoord(masses, bonds):
    """
    bonds: [ (index1, index2) ]
    """

    natom = len(masses)
    bondToAtoms = {}
    for ii in range(natom):
        bondToAtoms[ii] = []

    for ia, ja in bonds:
        bondToAtoms[ia].append(ja)
        bondToAtoms[ja].append(ia)

    terminal_atoms = _sort_atoms_by_mass([(ii, masses[ii])
                                          for ii in bondToAtoms
                                          if len(bondToAtoms[ii]) == 1],
                                         reverse=True)

    initial_atom = terminal_atoms[0][0]

    second_atom = bondToAtoms[initial_atom][0]

    list_atoms = []
    for iatom in bondToAtoms[second_atom]:
        if (iatom, masses[iatom]) in terminal_atoms:
            continue
        list_atoms.append((iatom, masses[iatom]))
    list_atoms.sort(key=sortFuncByMass, reverse=True)
    third_atom = list_atoms[0][0]

    root = [initial_atom, second_atom, third_atom]
    torsions = _find_torsion(root, masses, bondToAtoms)

    prior_atoms = [sorted([a1, a2]) for (a0, a1, a2, a3) in torsions]

    primary_torsion_indices = [prior_atoms.index(prior_atoms[n])
                               for n in range(len(prior_atoms))]

    unique_primary_torsion_indices = list(set(primary_torsion_indices))

    return torsions, primary_torsion_indices


def get_bat_from_xyz(xyz, zcoord, primary_torsion_indices):

    (a0, a1, a2, a3) = zcoord[0]

    p0 = xyz[a3]
    p1 = xyz[a2]
    p2 = xyz[a1]

    v01 = p1 - p0
    v21 = p1 - p2

    r01 = np.sqrt(np.sum(v01*v01))
    r12 = np.sqrt(np.sum(v21*v21))

    a012 = np.arccos(max(-1., min(1., np.sum(v01*v21) /
                                  np.sqrt(np.sum(v01*v01)*np.sum(v21*v21)))))
    e = v01 / r01
    phi = np.arctan2(e[1], e[0])  # Polar angle
    theta = np.arccos(e[2])  # Azimuthal angle
    # Rotation to the z axis
    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)
    Rz = np.array([[cp * ct, ct * sp, -st], [-sp, cp, 0],
                   [cp * st, sp * st, ct]])
    pos2 = Rz.dot(p2 - p1)

    # Angle about the rotation axis
    omega = np.arctan2(pos2[1], pos2[0])
    root_based = np.concatenate((p0, [phi, theta, omega, r01, r12, a012]))

    bonds = []
    angles = []
    torsions = []
    for a0, a1, a2, a3 in zcoord:
        r01 = _distance(xyz[a0], xyz[a1])
        a012 = _angle(xyz[a0], xyz[a1], xyz[a2])
        t0123 = _dihedral(xyz[a0], xyz[a1], xyz[a2], xyz[a3])

        bonds.append(r01)
        angles.append(a012)
        torsions.append(t0123)

    torsions = np.array(torsions)
    shift = torsions[primary_torsion_indices]
    unique_primary_torsion_indices = list(set(primary_torsion_indices))
    shift[unique_primary_torsion_indices] = 0.0

    torsions -= shift
    torsions = ((torsions + np.pi) % (2*np.pi)) - np.pi
    bat = np.concatenate((root_based, bonds, angles, torsions))

    return bat


def get_xyz_from_bat(bat, zcoord, primary_torsion_indices):

    # Split the bat vector inot more convenient variables
    origin = bat[:3]
    (phi, theta, omega) = bat[3:6]
    (r01, r12, a012) = bat[6:9]

    n_torsions = len(zcoord)

    bonds = bat[9:n_torsions+9]
    angles = bat[n_torsions+9:2*n_torsions+9]
    torsions = bat[2*n_torsions+9:]

    # convert improper to proper torsions
    shift = torsions[primary_torsion_indices]
    unique_torsion_indices = list(set(primary_torsion_indices))
    shift[unique_torsion_indices] = 0.0
    torsions += shift

    torsions = ((torsions+np.pi)) % (2.0*np.pi) - np.pi

    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([0.0, 0.0, r01])
    p2 = np.array([r12*np.sin(a012), 0, r01 - r12*np.cos(a012)])

    co = np.cos(omega)
    so = np.sin(omega)

    Romega = np.array([[co, -so, 0], [so, co, 0], [0, 0, 1]])
    p2 = Romega.dot(p2)

    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)
    # $R_Z(\phi) R_Y(\theta)$
    Re = np.array([[cp * ct, -sp, cp * st],
                   [ct * sp, cp, sp * st],
                   [-st, 0, ct]])
    p1 = Re.dot(p1)
    p2 = Re.dot(p2)
    # Translate the first three atoms by the origin
    p0 += origin
    p1 += origin
    p2 += origin

    xyz = np.zeros((n_torsions+3, 3), dtype=np.float32)
    (a0, a1, a2, a3) = zcoord[0]
    xyz[a3] = p0
    xyz[a2] = p1
    xyz[a1] = p2

    for ((a0, a1, a2, a3), r01, a012, t0123) \
            in zip(zcoord, bonds, angles, torsions):
        p1 = xyz[a1]
        p2 = xyz[a2]
        p3 = xyz[a3]

        sn_ang = np.sin(a012)
        cs_ang = np.cos(a012)
        sn_tor = np.sin(t0123)
        cs_tor = np.cos(t0123)

        v21 = p1 - p2
        v21 /= np.sqrt(np.sum(v21 * v21))
        v32 = p2 - p3
        v32 /= np.sqrt(np.sum(v32 * v32))

        vp = np.cross(v32, v21)
        cs = np.sum(v21 * v32)
        if abs(cs) > 1:
            print('cos ', cs)

        sn = np.sqrt(max(1.0 - cs * cs, 0.0000000001))
        vp = vp / sn
        vu = np.cross(vp, v21)

        xyz[a0] = p1 + \
            r01*(vu*sn_ang*cs_tor + vp*sn_ang*sn_tor - v21*cs_ang)

    return xyz


if __name__ == '__main__':

    import mdtraj

    traj = mdtraj.load('prmtopcrd/ligand.trans.inpcrd',
                       top='prmtopcrd/ligand.prmtop')

    xyz = traj.xyz[0]*10.0

    print(xyz)

    natom = traj.topology._numAtoms

    masses = [a.element.mass for a in traj.topology.atoms]
    bonds = [(b.atom1.index, b.atom2.index)
             for b in traj.topology.bonds]

    zcoord, primary_torsion_indices = get_zcoord(masses, bonds)

    print(len(zcoord))
    print(primary_torsion_indices, len(primary_torsion_indices))
    unique_primary_torsion_indices = list(set(primary_torsion_indices))
    print(unique_primary_torsion_indices)

    bat = get_bat_from_xyz(
        xyz, zcoord, primary_torsion_indices)

    xyz1 = get_xyz_from_bat(bat, zcoord, primary_torsion_indices)
    print(xyz1)
