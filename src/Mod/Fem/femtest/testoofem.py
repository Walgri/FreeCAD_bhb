# ***************************************************************************
# *   Copyright (c) 2019 Bernd Hahnebach <bernd@bimstatik.org>              *
# *                                                                         *
# *   This file is part of the FreeCAD CAx development system.              *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   FreeCAD is distributed in the hope that it will be useful,            *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with FreeCAD; if not, write to the Free Software        *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************


import FreeCAD
import unittest
from . import utilstest as testtools
from .utilstest import fcc_print

from os.path import join


class TestOofem(unittest.TestCase):
    fcc_print('import TestOofem')

    def setUp(self):
        self.doc_name = "TestOofem"
        try:
            FreeCAD.setActiveDocument(self.doc_name)
        except:
            FreeCAD.newDocument(self.doc_name)
        finally:
            FreeCAD.setActiveDocument(self.doc_name)
        self.active_doc = FreeCAD.ActiveDocument
        self.mesh_name = 'Mesh'
        self.temp_dir = testtools.get_fem_test_tmp_dir()
        self.test_file_dir = join(testtools.get_fem_test_home_dir(), 'oofem')
        fcc_print('\nDocument: {} setup for testing'.format(self.active_doc.Name))

    def test_planestress2d(self):
        # test oofem input file writing for oofem manual example PlaneStress2D.in
        from femexamples.oofem_planestress2d import setup_planestress2d as planestress2d
        planestress2d(self.active_doc, 'oofem')
        self.write_inputfile('OOFEM_PlaneStress2D')

    def test_cantileverccxnodeload(self):
        # test oofem input file writing for CalculiX cantilever loaded by 9 MN as node loads
        from femexamples.ccx_cantilever_std import setup_cantilevernodeload as cantinoadload
        cantinoadload(self.active_doc, 'oofem')
        self.write_inputfile('OOFEM_CantileverCcxNodeLoad')

    def test_cantileverccxprescribeddisplacement(self):
        # test oofem input file writing for CalculiX cantilever loaded by 9 MN as node loads
        from femexamples.ccx_cantilever_std import setup_cantileverprescribeddisplacement as cantipredisp
        cantipredisp(self.active_doc, 'oofem')
        self.write_inputfile('OOFEM_CantileverCcxPrescribedDisplacement')

    def write_inputfile(self, static_base_name):

        # recompute
        self.active_doc.recompute()

        # fcc_print(self.active_doc.Objects)
        # fcc_print([obj.Name for obj in self.active_doc.Objects])

        # save file, needs before input file writing because of besides dir
        solverframework_analysis_dir = testtools.get_unit_test_tmp_dir(testtools.get_fem_test_tmp_dir(), 'FEM_oofem')
        save_fc_file = join(solverframework_analysis_dir, static_base_name) + '.FCStd'
        fcc_print('Save FreeCAD file for {} analysis to {}'.format(static_base_name, save_fc_file))
        self.active_doc.saveAs(save_fc_file)

        # solver
        solver = self.active_doc.SolverOOFEM

        # get analysis workig dir
        from femsolver.run import _getBesideDir as getpath
        working_dir = getpath(solver)

        fcc_print('Checking FEM OOFEM solver ...')
        fcc_print('machine_oofem')

        # write input file
        from femsolver.run import PREPARE
        # set up the machine with testmode True
        machine_oofem = solver.Proxy.createMachine(solver, working_dir, True)
        machine_oofem.target = PREPARE
        machine_oofem.start()
        machine_oofem.join()  # wait for the machine to finish.

        # compare input files
        inpfile_given = join(testtools.get_fem_test_home_dir(), 'oofem', (static_base_name + '.in'))
        inpfile_totest = join(working_dir, (self.mesh_name + '.in'))
        fcc_print('Comparing {}  to  {}'.format(inpfile_given, inpfile_totest))
        ret = testtools.compare_inp_files(inpfile_given, inpfile_totest)
        self.assertFalse(ret, "OOFEM write_inp_file for {0} test failed.\n{1}".format(static_base_name, ret))

    def tearDown(self):
        FreeCAD.closeDocument(self.doc_name)
        fcc_print('Document closed')
        pass


'''
make -j4  &&  ./bin/FreeCAD --run-test "femtest.testoofem.TestOofem"
make -j4  &&  ./bin/FreeCAD --run-test "femtest.testoofem.TestOofem.test_planestress2d"
make -j4  &&  ./bin/FreeCAD --run-test "femtest.testoofem.TestOofem.test_cantileverccxnodeload"
make -j4  &&  ./bin/FreeCAD --run-test "femtest.testoofem.TestOofem.test_cantileverccxprescribeddisplacement"

'''
