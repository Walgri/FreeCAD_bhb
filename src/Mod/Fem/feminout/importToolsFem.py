# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2017 - Bernd Hahnebach <bernd@bimstatik.org>            *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

__title__ = "FreeCAD FEM import tools"
__author__ = "Bernd Hahnebach"
__url__ = "http://www.freecadweb.org"

## @package importToolsFem
#  \ingroup FEM
#  \brief FreeCAD FEM import tools

import numpy as np

import FreeCAD
import femmesh.meshtools


def get_FemMeshObjectMeshGroups(fem_mesh_obj):
    """
        Get mesh groups from mesh. This also throws no exception if there
        is no Groups property at all (e.g. Netgen meshes).
    """
    fem_mesh = fem_mesh_obj.FemMesh
    try:
        gmshgroups = fem_mesh.Groups
    except:
        gmshgroups = ()

    return gmshgroups


def get_FemMeshObjectOrder(fem_mesh_obj):
    """
        Gets element order. Element order counting based on number of nodes on
        edges. Edge with 2 nodes -> linear elements, Edge with 3 nodes ->
        quadratic elements, and so on. No edges in mesh -> not determined.
        (Is this possible? Seems to be a very degenerate case.)
        If there are edges with different number of nodes appearing, return
        list of orders.
    """
    presumable_order = None

    edges = fem_mesh_obj.FemMesh.Edges

    if edges != ():
        edges_length_set = list({len(fem_mesh_obj.FemMesh.getElementNodes(e)) for e in edges})
        # only need set to eliminate double entries

        if len(edges_length_set) == 1:
            presumable_order = edges_length_set[0] - 1
        else:
            presumable_order = [el - 1 for el in edges_length_set]
    else:
        print("Found no edges in mesh: Element order determination does not work without them.")

    return presumable_order


def get_FemMeshObjectDimension(fem_mesh_obj):
    """ Count all entities in an abstract sense, to distinguish which dimension the mesh is
        (i.e. linemesh, facemesh, volumemesh)
    """
    dim = None

    if fem_mesh_obj.FemMesh.Nodes != ():
        dim = 0
    if fem_mesh_obj.FemMesh.Edges != ():
        dim = 1
    if fem_mesh_obj.FemMesh.Faces != ():
        dim = 2
    if fem_mesh_obj.FemMesh.Volumes != ():
        dim = 3

    return dim


def get_FemMeshObjectElementTypes(fem_mesh_obj, remove_zero_element_entries=True):
    """
        Spit out all elements in the mesh with their appropriate dimension.
    """
    FreeCAD_element_names_dims = {
        "Node": 0, "Edge": 1, "Hexa": 3, "Polygon": 2, "Polyhedron": 3,
        "Prism": 3, "Pyramid": 3, "Quadrangle": 2, "Tetra": 3, "Triangle": 2}

    eval_dict = locals()  # to access local variables from eval
    elements_list_with_zero = [(eval("fem_mesh_obj.FemMesh." + s + "Count", eval_dict), s, d) for (s, d) in FreeCAD_element_names_dims.items()]
    # ugly but necessary
    if remove_zero_element_entries:
        elements_list = [(num, s, d) for (num, s, d) in elements_list_with_zero if num > 0]
    else:
        elements_list = elements_list_with_zero

    return elements_list


def get_MaxDimElementFromList(elem_list):
    """
        Gets element with the maximal dimension in the mesh to determine cells.
    """
    elem_list.sort(key=lambda t: t[2])
    return elem_list[-1]


def make_femmesh(mesh_data):
    ''' makes an FreeCAD FEM Mesh object from FEM Mesh data
    '''
    import Fem
    mesh = Fem.FemMesh()
    m = mesh_data
    if ('Nodes' in m) and (len(m['Nodes']) > 0):
        FreeCAD.Console.PrintLog("Found: nodes\n")
        if (
            ('Seg2Elem' in m)
            or ('Seg3Elem' in m)
            or ('Tria3Elem' in m)
            or ('Tria6Elem' in m)
            or ('Quad4Elem' in m)
            or ('Quad8Elem' in m)
            or ('Tetra4Elem' in m)
            or ('Tetra10Elem' in m)
            or ('Penta6Elem' in m)
            or ('Penta15Elem' in m)
            or ('Hexa8Elem' in m)
            or ('Hexa20Elem' in m)
        ):

            nds = m['Nodes']
            FreeCAD.Console.PrintLog("Found: elements\n")
            for i in nds:
                n = nds[i]
                mesh.addNode(n[0], n[1], n[2], i)
            elms_hexa8 = m['Hexa8Elem']
            for i in elms_hexa8:
                e = elms_hexa8[i]
                mesh.addVolume([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]], i)
            elms_penta6 = m['Penta6Elem']
            for i in elms_penta6:
                e = elms_penta6[i]
                mesh.addVolume([e[0], e[1], e[2], e[3], e[4], e[5]], i)
            elms_tetra4 = m['Tetra4Elem']
            for i in elms_tetra4:
                e = elms_tetra4[i]
                mesh.addVolume([e[0], e[1], e[2], e[3]], i)
            elms_tetra10 = m['Tetra10Elem']
            for i in elms_tetra10:
                e = elms_tetra10[i]
                mesh.addVolume([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9]], i)
            elms_penta15 = m['Penta15Elem']
            for i in elms_penta15:
                e = elms_penta15[i]
                mesh.addVolume([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9],
                                e[10], e[11], e[12], e[13], e[14]], i)
            elms_hexa20 = m['Hexa20Elem']
            for i in elms_hexa20:
                e = elms_hexa20[i]
                mesh.addVolume([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9],
                                e[10], e[11], e[12], e[13], e[14], e[15], e[16], e[17], e[18], e[19]], i)
            elms_tria3 = m['Tria3Elem']
            for i in elms_tria3:
                e = elms_tria3[i]
                mesh.addFace([e[0], e[1], e[2]], i)
            elms_tria6 = m['Tria6Elem']
            for i in elms_tria6:
                e = elms_tria6[i]
                mesh.addFace([e[0], e[1], e[2], e[3], e[4], e[5]], i)
            elms_quad4 = m['Quad4Elem']
            for i in elms_quad4:
                e = elms_quad4[i]
                mesh.addFace([e[0], e[1], e[2], e[3]], i)
            elms_quad8 = m['Quad8Elem']
            for i in elms_quad8:
                e = elms_quad8[i]
                mesh.addFace([e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]], i)
            elms_seg2 = m['Seg2Elem']
            for i in elms_seg2:
                e = elms_seg2[i]
                mesh.addEdge([e[0], e[1]], i)
            elms_seg3 = m['Seg3Elem']
            for i in elms_seg3:
                e = elms_seg3[i]
                mesh.addEdge([e[0], e[1], e[2]], i)
            FreeCAD.Console.PrintLog("imported mesh: {} nodes, {} HEXA8, {} PENTA6, {} TETRA4, {} TETRA10, {} PENTA15".format(
                len(nds), len(elms_hexa8), len(elms_penta6), len(elms_tetra4), len(elms_tetra10), len(elms_penta15)
            ))
            FreeCAD.Console.PrintLog("imported mesh: {} HEXA20, {} TRIA3, {} TRIA6, {} QUAD4, {} QUAD8, {} SEG2, {} SEG3".format(
                len(elms_hexa20), len(elms_tria3), len(elms_tria6), len(elms_quad4), len(elms_quad8), len(elms_seg2), len(elms_seg3)
            ))
        else:
            FreeCAD.Console.PrintError("No Elements found!\n")
    else:
        FreeCAD.Console.PrintError("No Nodes found!\n")
    return mesh


def fill_femresult_mechanical(res_obj, result_set):
    ''' fills a FreeCAD FEM mechanical result object with result data
    '''
    if 'time' in result_set:
        step_time = result_set['time']
        step_time = round(step_time, 2)

    # if disp exists, fill res_obj.NodeNumbers and res_obj.DisplacementVectors as well as stress and strain
    if 'disp' in result_set:
        disp = result_set['disp']
        displacement = []
        for k, v in disp.items():
            displacement.append(v)

        x_max, y_max, z_max = map(max, zip(*displacement))
        if eigenmode_number > 0:
            span = get_span(res_obj.Mesh.FemMesh.Nodes.items())
            max_disp = max(x_max, y_max, z_max)
            # Allow for max displacement to be 0.1% of the span
            # FIXME - add to Preferences
            max_allowed_disp = 0.001 * span
            scale = max_allowed_disp / max_disp
            res_obj.EigenmodeScaleFactor = scale
            scale = 1.0

        res_obj.DisplacementVectors = list(map((lambda x: x * scale), disp.values()))
        res_obj.NodeNumbers = list(disp.keys())

        # fill res_obj.NodeStressXX etc if they exist in result_set
        # list values are just added
        # Should we gone check if the key in stress and strain dict is the same as the number in NodeNumbers?
        if 'stress' in result_set:
            stress = result_set['stress']
            Sxx = []
            Syy = []
            Szz = []
            Sxy = []
            Syz = []
            Szx = []
            for i, values_S in enumerate(stress.values()):  # values_S .. stresstuple .. (Sxx, Syy, Szz, Sxy, Syz, Szx)
                Sxx.append(values_S[0])
                Syy.append(values_S[1])
                Szz.append(values_S[2])
                Sxy.append(values_S[3])
                Syz.append(values_S[4])
                Szx.append(values_S[5])
            res_obj.NodeStressXX = Sxx
            res_obj.NodeStressYY = Syy
            res_obj.NodeStressZZ = Szz
            res_obj.NodeStressXY = Sxy
            res_obj.NodeStressYZ = Syz
            res_obj.NodeStressZX = Szx

        # fill res_obj.NodeStrainXX etc if they exist in result_set
        if 'strain' in result_set:
            strain = result_set['strain']
            Exx = []
            Eyy = []
            Ezz = []
            Exy = []
            Eyz = []
            Ezx = []
            for i, values_E in enumerate(strain.values()):  # values_E .. straintuple .. (Exx, Eyy, Ezz, Exy, Eyz, Ezx)
                Exx.append(values_E[0])
                Eyy.append(values_E[1])
                Ezz.append(values_E[2])
                Exy.append(values_E[3])
                Eyz.append(values_E[4])
                Ezx.append(values_E[5])
            res_obj.NodeStrainXX = Exx
            res_obj.NodeStrainYY = Eyy
            res_obj.NodeStrainZZ = Ezz
            res_obj.NodeStrainXY = Exy
            res_obj.NodeStrainYZ = Eyz
            res_obj.NodeStrainZX = Ezx

        # calculate principal and max Shear and fill them in res_obj
        if 'stress' in result_set:
            stress = result_set['stress']
            nsr = len(stress)
            print("nsr: {}".format(nsr))

            if nsr > 0:
                prinstress1 = []
                prinstress2 = []
                prinstress3 = []
                shearstress = []
                ps1v = []
                ps2v = []
                ps3v = []
                ic = np.zeros(nsr)
                #
                # HarryvL: addtional arrays to hold reinforcement ratios
                # and mohr coulomb stress
                #
                rhx = []
                rhy = []
                rhz = []
                moc = []
                #
                # HarryvL: determine concrete / non-concrete nodes
                #
                result_mesh = res_obj.Mesh.FemMesh
                for obj in FreeCAD.ActiveDocument.Objects:
                    if obj.isDerivedFrom('App::MaterialObjectPython'):
                        if obj.Material.get('Name') == "Concrete":
                            print("CONCRETE")
                            if obj.References == []:
                                for iic in range(nsr):
                                    if ic[iic] == 0:
                                        ic[iic] = 1
                            else:
                                for ref in obj.References:
                                    concrete_nodes = femmesh.meshtools.get_femnodes_by_refshape(result_mesh, ref)
                                    for cn in concrete_nodes:
                                        ic[cn - 1] = 1
                        else:
                            print("NOT CONCRETE")
                            if obj.References == []:
                                for iic in range(nsr):
                                    if ic[iic] == 0:
                                        ic[iic] = 2
                            else:
                                for ref in obj.References:
                                    non_concrete_nodes = femmesh.meshtools.get_femnodes_by_refshape(result_mesh, ref)
                                    for ncn in non_concrete_nodes:
                                        ic[ncn - 1] = 2

                for isv in range(nsr):

                    i = list(stress.values())[isv]

                    rhox = 0.
                    rhoy = 0.
                    rhoz = 0.
                    mc = 0.
                    scxx = i[0]
                    scyy = i[1]
                    sczz = i[2]

                    if ic[isv] == 1:
                        #
                        # HarryvL: for concrete scxx etc. are affected by
                        # reinforcement (see calculate_rho(i)). for all other
                        # materials scxx etc. are the original stresses
                        #
                        rhox, rhoy, rhoz, scxx, scyy, sczz = calculate_rho(i)

                    prin1, prin2, prin3, shear, psv =\
                        calculate_principal_stress(i, scxx, scyy, sczz)

                    prinstress1.append(prin1)
                    prinstress2.append(prin2)
                    prinstress3.append(prin3)
                    shearstress.append(shear)
                    ps1v.append(psv[0])
                    ps2v.append(psv[1])
                    ps3v.append(psv[2])

                    #
                    # reinforcement ratios and mohr coulomb criterion
                    #
                    rhx.append(rhox)
                    rhy.append(rhoy)
                    rhz.append(rhoz)
                    if ic[isv] == 1:
                        mc = calculate_mohr_coulomb(prin1, prin3)
                    moc.append(mc)

                if eigenmode_number > 0:
                    res_obj.PrincipalMax = list(map((lambda x: x * scale),
                                                    prinstress1))
                    res_obj.PrincipalMed = list(map((lambda x: x * scale),
                                                    prinstress2))
                    res_obj.PrincipalMin = list(map((lambda x: x * scale),
                                                    prinstress3))
                    res_obj.MaxShear = list(map((lambda x: x * scale),
                                                shearstress))

                    #
                    # HarryvL: addtional concrete and principal stress
                    # results for use in _ViewProviderFemResultMechanical
                    #
                    res_obj.ReinforcementRatio_x = list(map((lambda x:
                                                             x * scale), rhx))
                    res_obj.ReinforcementRatio_y = list(map((lambda x:
                                                             x * scale), rhy))
                    res_obj.ReinforcementRatio_z = list(map((lambda x:
                                                             x * scale), rhz))
                    res_obj.MohrCoulomb = list(map((lambda x: x * scale), moc))

                    res_obj.PS1Vector = list(map((lambda x: x * scale), ps1v))
                    res_obj.PS2Vector = list(map((lambda x: x * scale), ps2v))
                    res_obj.PS3Vector = list(map((lambda x: x * scale), ps3v))

                    res_obj.Eigenmode = eigenmode_number
                else:
                    res_obj.PrincipalMax = prinstress1
                    res_obj.PrincipalMed = prinstress2
                    res_obj.PrincipalMin = prinstress3
                    res_obj.MaxShear = shearstress
                    #
                    # HarryvL: addtional concrete and principal stress plot
                    # results for use in _ViewProviderFemResultMechanical
                    #
                    res_obj.ReinforcementRatio_x = rhx
                    res_obj.ReinforcementRatio_y = rhy
                    res_obj.ReinforcementRatio_z = rhz
                    res_obj.MohrCoulomb = moc

                    res_obj.PS1Vector = ps1v
                    res_obj.PS2Vector = ps2v
                    res_obj.PS3Vector = ps3v

        # fill Equivalent Plastic strain if they exist
        if 'peeq' in result_set:
            Peeq = result_set['peeq']
            if len(Peeq) > 0:
                if len(Peeq.values()) != len(disp.values()):
                    Pe = []
                    Pe_extra_nodes = Peeq.values()
                    nodes = len(disp.values())
                    for i in range(nodes):
                        Pe_value = Pe_extra_nodes[i]
                        Pe.append(Pe_value)
                    res_obj.Peeq = Pe
                else:
                    res_obj.Peeq = Peeq.values()

    # fill res_obj.Temperature if they exist
    # TODO, check if it is possible to have Temperature without disp, we would need to set NodeNumbers than
    if 'temp' in result_set:
        Temperature = result_set['temp']
        if len(Temperature) > 0:
            if len(Temperature.values()) != len(disp.values()):
                Temp = []
                Temp_extra_nodes = Temperature.values()
                nodes = len(disp.values())
                for i in range(nodes):
                    Temp_value = Temp_extra_nodes[i]
                    Temp.append(Temp_value)
                res_obj.Temperature = list(map((lambda x: x), Temp))
            else:
                res_obj.Temperature = list(map((lambda x: x), Temperature.values()))
            res_obj.Time = step_time

    # fill res_obj.MassFlow
    if 'mflow' in result_set:
        MassFlow = result_set['mflow']
        if len(MassFlow) > 0:
            res_obj.MassFlowRate = list(map((lambda x: x), MassFlow.values()))
            res_obj.Time = step_time
            res_obj.NodeNumbers = list(MassFlow.keys())  # disp does not exist, res_obj.NodeNumbers needs to be set

    # fill res_obj.NetworkPressure, disp does not exist, see MassFlow
    if 'npressure' in result_set:
        NetworkPressure = result_set['npressure']
        if len(NetworkPressure) > 0:
            res_obj.NetworkPressure = list(map((lambda x: x), NetworkPressure.values()))
            res_obj.Time = step_time

    return res_obj


# helper
def get_span(node_items):
    positions = []  # list of node vectors
    for k, v in node_items:
        positions.append(v)
    p_x_max, p_y_max, p_z_max = map(max, zip(*positions))
    p_x_min, p_y_min, p_z_min = map(min, zip(*positions))
    x_span = abs(p_x_max - p_x_min)
    y_span = abs(p_y_max - p_y_min)
    z_span = abs(p_z_max - p_z_min)
    span = max(x_span, y_span, z_span)
    return span


def calculate_principal_stress(i, scxx, scyy, sczz):
    #
    #   HarryvL - calculate principal stress vectors and values
    #           - for concrete stresses use scxx, scyy, sczz on the diagonal
    #             of the stress tensor
    #           - for total stresses use i[0], i[1], i[2] on the diagonal of
    #             the stress tensor
    #           - TODO: option to use concrete or total stresses by user
    #
    #

    sigma = np.array([[i[0], i[3], i[5]],
                      [i[3], i[1], i[4]],
                      [i[5], i[4], i[3]]])

    eigenvalues, eigenvectors = np.linalg.eig(sigma)

    #
    #   HarryvL: suppress complex eigenvalue and vectors that may occur for
    #   near-zero (numerical noise) stress fields
    #

    eigenvalues = eigenvalues.real
    eigenvectors = eigenvectors.real

    eigenvectors[:, 0] = eigenvalues[0] * eigenvectors[:, 0]
    eigenvectors[:, 1] = eigenvalues[1] * eigenvectors[:, 1]
    eigenvectors[:, 2] = eigenvalues[2] * eigenvectors[:, 2]

    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    maxshear = (eigenvalues[0] - eigenvalues[2]) / 2.0

    return (eigenvalues[0], eigenvalues[1], eigenvalues[2], maxshear,
            tuple([tuple(row) for row in eigenvectors.T]))


def calculate_rho(i):

    #
    #   HarryvL - Calculation of Reinforcement Ratios and
    #   Concrete Stresses according to http://heronjournal.nl/53-4/3.pdf
    #           - See post:
    #             https://forum.freecadweb.org/viewtopic.php?f=18&t=28821
    #           - TODO: the following material parameters are hard-coded
    #             and should be entered in material dialog
    #                   fy: factored yield strength of reinforcement bars
    #                   r0: optional value for minimum reinforcement ratio
    #                       (rx, y, rz are reduced accordingly) - default 0.0
    #
    fy = 315.
    r0 = 0.0

    rmin = 1.0e9
    eqmin = 14

    sxx = i[0]
    syy = i[1]
    szz = i[2]
    sxy = i[3]
    syz = i[4]
    sxz = i[5]

    rhox = np.zeros(15)
    rhoy = np.zeros(15)
    rhoz = np.zeros(15)

    #    i1=sxx+syy+szz NOT USED
    #    i2=sxx*syy+syy*szz+szz*sxx-sxy**2-sxz**2-syz**2 NOT USED
    i3 = (sxx * syy * szz + 2 * sxy * sxz * syz - sxx * syz**2
          - syy * sxz**2 - szz * sxy**2)

    #    Solution (5)
    d = (sxx * syy - sxy**2)
    if d != 0.:
        rhoz[0] = i3 / d / fy

    #    Solution (6)
    d = (sxx * szz - sxz**2)
    if d != 0.:
        rhoy[1] = i3 / d / fy

    #    Solution (7)
    d = (syy * szz - syz**2)
    if d != 0.:
        rhox[2] = i3 / d / fy

    #    Solution (9)
    if sxx != 0.:
        fc = sxz * sxy / sxx - syz
        fxy = sxy**2 / sxx
        fxz = sxz**2 / sxx

        #    Solution (9+)
        rhoy[3] = syy - fxy + fc
        rhoy[3] /= fy
        rhoz[3] = szz - fxz + fc
        rhoz[3] /= fy

        #    Solution (9-)
        rhoy[4] = syy - fxy - fc
        rhoy[4] /= fy
        rhoz[4] = szz - fxz - fc
        rhoz[4] /= fy

    #   Solution (10)
    if syy != 0.:
        fc = syz * sxy / syy - sxz
        fxy = sxy**2 / syy
        fyz = syz**2 / syy

        # Solution (10+)
        rhox[5] = sxx - fxy + fc
        rhox[5] /= fy
        rhoz[5] = szz - fyz + fc
        rhoz[5] /= fy

        # Solution (10-)vm
        rhox[6] = sxx - fxy - fc

        rhox[6] /= fy
        rhoz[6] = szz - fyz - fc
        rhoz[6] /= fy

    # Solution (11)
    if szz != 0.:
        fc = sxz * syz / szz - sxy
        fxz = sxz**2 / szz
        fyz = syz**2 / szz

        # Solution (11+)
        rhox[7] = sxx - fxz + fc
        rhox[7] /= fy
        rhoy[7] = syy - fyz + fc
        rhoy[7] /= fy

        # Solution (11-)
        rhox[8] = sxx - fxz - fc
        rhox[8] /= fy
        rhoy[8] = syy - fyz - fc
        rhoy[8] /= fy

    # Solution (13)
    rhox[9] = (sxx + sxy + sxz) / fy
    rhoy[9] = (syy + sxy + syz) / fy
    rhoz[9] = (szz + sxz + syz) / fy

    # Solution (14)
    rhox[10] = (sxx + sxy - sxz) / fy
    rhoy[10] = (syy + sxy - syz) / fy
    rhoz[10] = (szz - sxz - syz) / fy

    # Solution (15)
    rhox[11] = (sxx - sxy - sxz) / fy
    rhoy[11] = (syy - sxy + syz) / fy
    rhoz[11] = (szz - sxz + syz) / fy

    # Solution (16)
    rhox[12] = (sxx - sxy + sxz) / fy
    rhoy[12] = (syy - sxy - syz) / fy
    rhoz[12] = (szz + sxz - syz) / fy

    # Solution (17)
    if syz != 0.:
        rhox[13] = (sxx - sxy * sxz / syz) / fy
    if sxz != 0.:
        rhoy[13] = (syy - sxy * syz / sxz) / fy
    if sxy != 0.:
        rhoz[13] = (szz - sxz * syz / sxy) / fy

    for ir in range(0, rhox.size):

        if rhox[ir] >= -1.e-10 and rhoy[ir] >= -1.e-10 and rhoz[ir] > -1.e-10:

            # Concrete Stresses
            scxx = sxx - rhox[ir] * fy
            scyy = syy - rhoy[ir] * fy
            sczz = szz - rhoz[ir] * fy
            ic1 = (scxx + scyy + sczz)
            ic2 = (scxx * scyy + scyy * sczz + sczz * scxx - sxy**2
                   - sxz**2 - syz**2)
            ic3 = (scxx * scyy * sczz + 2 * sxy * sxz * syz - scxx * syz**2
                   - scyy * sxz**2 - sczz * sxy**2)

            if ic1 <= 1.e-6 and ic2 >= -1.e-6 and ic3 <= 1.0e-6:

                rsum = rhox[ir] + rhoy[ir] + rhoz[ir]

                if rsum < rmin and rsum > 0.:
                    rmin = rsum
                    eqmin = ir

    rx = max(rhox[eqmin] - r0, 0.0)
    ry = max(rhoy[eqmin] - r0, 0.0)
    rz = max(rhoz[eqmin] - r0, 0.0)

    scxx = sxx - rx * fy
    scyy = syy - ry * fy
    sczz = szz - rz * fy

    return rhox[eqmin], rhoy[eqmin], rhoz[eqmin], scxx, scyy, sczz


def calculate_mohr_coulomb(prin1, prin3):
    #
    #   HarryvL - Calculation of Mohr Coulomb yield criterion to judge
    #             concrete curshing and shear failure
    #           - TODO: the following material parameters are hard-coded
    #             and should be entered in material dialog
    #                   phi: angle of internal friction for
    #                        concrete - default 30 degrees
    #                   fck: factored concrete cube compressive
    #                        stength - default 0.75*0.6*35.0 = 15.75 MPa
    #

    phi = np.pi / 6.
    fck = 15.75
    coh = fck * (1 - np.sin(phi)) / 2 / np.cos(phi)

    mc_stress = ((prin1 - prin3) + (prin1 + prin3) * np.sin(phi)
                 - 2. * coh * np.cos(phi))

    if mc_stress < 0.:
        mc_stress = 0.

    return mc_stress
