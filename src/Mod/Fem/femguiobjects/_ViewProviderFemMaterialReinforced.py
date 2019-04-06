# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2019 Bernd Hahnebach <bernd@bimstatik.org>              *
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

__title__ = "FreeCAD FEM material reinforced ViewProvider for the document object"
__author__ = "Bernd Hahnebach"
__url__ = "http://www.freecadweb.org"

## @package ViewProviderFemMaterialReinforced
#  \ingroup FEM
#  \brief FreeCAD FEM _ViewProviderFemMaterialReinforced

import FreeCAD
import FreeCADGui
import FemGui  # needed to display the icons in TreeView
False if False else FemGui.__name__  # flake8, dummy FemGui usage, returns 'FemGui'

# task panel
# from . import FemSelectionWidgets
from PySide import QtCore
from PySide import QtGui
import sys
import MaterialEditor


class _ViewProviderFemMaterialReinforced:
    "A View Provider for the FemMaterialReinfocement object"
    def __init__(self, vobj):
        vobj.Proxy = self

    def getIcon(self):
        return ":/icons/fem-material-reinforced.svg"

    def attach(self, vobj):
        from pivy import coin
        self.ViewObject = vobj
        self.Object = vobj.Object
        self.standard = coin.SoGroup()
        vobj.addDisplayMode(self.standard, "Default")

    def getDisplayModes(self, obj):
        return ["Default"]

    def updateData(self, obj, prop):
        return

    def onChanged(self, vobj, prop):
        return

    def setEdit(self, vobj, mode=0):
        # hide all meshes
        for o in FreeCAD.ActiveDocument.Objects:
            if o.isDerivedFrom("Fem::FemMeshObject"):
                o.ViewObject.hide()
        # hide all meshes
        for o in FreeCAD.ActiveDocument.Objects:
            if o.isDerivedFrom("Fem::FemMeshObject"):
                o.ViewObject.hide()
        # show task panel
        taskd = _TaskPanelFemMaterialReinforced(self.Object)
        taskd.obj = vobj.Object
        FreeCADGui.Control.showDialog(taskd)
        return False

    def unsetEdit(self, vobj, mode=0):
        FreeCADGui.Control.closeDialog()
        return True

    def doubleClicked(self, vobj):
        guidoc = FreeCADGui.getDocument(vobj.Object.Document)
        # check if another VP is in edit mode
        # https://forum.freecadweb.org/viewtopic.php?t=13077#p104702
        if not guidoc.getInEdit():
            guidoc.setEdit(vobj.Object.Name)
        else:
            from PySide.QtGui import QMessageBox
            message = (
                'Active Task Dialog found! '
                'Please close this one before opening  a new one!'
            )
            QMessageBox.critical(None, "Error in tree view", message)
            FreeCAD.Console.PrintError(message + '\n')
        return True

    def __getstate__(self):
        return None

    def __setstate__(self, state):
        return None


class _TaskPanelFemMaterialReinforced:
    '''The editmode TaskPanel for FemMaterialReinforced objects'''

    if sys.version_info.major >= 3:
        unicode = str

    def __init__(self, obj):

        self.obj = obj
        self.output_obj_mat_param()

        # init matrix and reinforcement material
        self.matrix = self.obj.Material
        self.card_path_m = ''
        self.has_transient_mat_m = False
        self.reinforcement = self.obj.Reinforcement
        self.card_path_r = ''
        self.has_transient_mat_r = False
        # mat_card is the FCMat file
        # card_name is the file name of the mat_card
        # card_path is the whole file path of the mat_card
        # material_name is the value of the key name in FreeCAD material dictionary
        # they might not match because of special letters in the material_name which are
        # changed in the card_name to english standard characters

        # init for collecting all mat data and icons
        self.materials = {}  # { card_path : FreeCAD material dict }
        self.icons = {}  # { card_path : icon_path }

        # parameter widget
        self.parameterWidget = FreeCADGui.PySideUic.loadUi(
            FreeCAD.getHomePath() + "Mod/Fem/Resources/ui/MaterialReinforcement.ui"
        )

        # globals
        QtCore.QObject.connect(
            self.parameterWidget.cb_materials_m,
            QtCore.SIGNAL("activated(int)"),
            self.choose_material_m
        )
        QtCore.QObject.connect(
            self.parameterWidget.pb_edit_m,
            QtCore.SIGNAL("clicked()"),
            self.edit_material_m
        )
        QtCore.QObject.connect(
            self.parameterWidget.cb_materials_r,
            QtCore.SIGNAL("activated(int)"),
            self.choose_material_r
        )
        QtCore.QObject.connect(
            self.parameterWidget.pb_edit_r,
            QtCore.SIGNAL("clicked()"),
            self.edit_material_r
        )

        # get all available materials
        self.import_materials()
        # fill the matrix and reinforcement material comboboxes with material cards
        self.add_cards_to_combo_boxes()

        # search for exact the mat_card_m and mat_card_r in all known cards
        # choose the current matrix material
        self.card_path_m = self.get_material_card(self.matrix)
        print('card_path: ' + self.card_path_m)
        if not self.card_path_m:
            # we have not found our material in self.materials dict :-(
            # we're going to add a user-defined temporary material: a document material
            FreeCAD.Console.PrintMessage(
                "Previously used material card cannot be found in material directories. "
                "Add document material.\n"
            )
            self.card_path_m = '_Document_Matrix_Materialx'
            self.materials[self.card_path_m] = self.matrix
            self.parameterWidget.cb_materials_m.addItem(
                QtGui.QIcon(":/icons/help-browser.svg"),
                self.card_path_m,
                self.card_path_m
            )
            index = self.parameterWidget.cb_materials_m.findData(self.card_path_m)
            # print(index)
            # fill input fields and set the current material in the cb widget
            self.choose_material_m(index)
        else:
            # we found our exact material in self.materials dict :-)
            FreeCAD.Console.PrintMessage(
                "Previously used material card was found in material directories. "
                "We will use this material.\n"
            )
            index = self.parameterWidget.cb_materials_m.findData(self.card_path_m)
            # set the current material in the cb widget
            self.choose_material_m(index)

        # choose the current reinforcement material
        self.card_path_r = self.get_material_card(self.reinforcement)
        print('card_path: ' + self.card_path_r)
        if not self.card_path_r:
            # we have not found our material in self.materials dict :-(
            # we're going to add a user-defined temporary material: a document material
            FreeCAD.Console.PrintMessage(
                "Previously used material card cannot be found in material directories. "
                "Add document material.\n"
            )
            self.card_path_r = '_Document_Reinforcement_Material'
            self.materials[self.card_path_r] = self.reinforcement
            self.parameterWidget.cb_materials_r.addItem(
                QtGui.QIcon(":/icons/help-browser.svg"),
                self.card_path_r,
                self.card_path_r
            )
            index = self.parameterWidget.cb_materials_r.findData(self.card_path_r)
            # set the current material in the cb widget
            self.choose_material_r(index)
        else:
            # we found our exact material in self.materials dict :-)
            FreeCAD.Console.PrintMessage(
                "Previously used material card was found in material directories. "
                "We will use this material.\n"
            )
            index = self.parameterWidget.cb_materials_r.findData(self.card_path_r)
            # print(index)
            # fill input fields and set the current material in the cb widget
            self.choose_material_r(index)

        # set up the form
        self.form = self.parameterWidget

    # leave task panel ***************************************************************************
    def accept(self):
        self.obj.Material = self.matrix
        self.obj.Reinforcement = self.reinforcement
        self.recompute_and_set_back_all()
        return True

    def reject(self):
        self.recompute_and_set_back_all()
        return True

    def recompute_and_set_back_all(self):
        guidoc = FreeCADGui.getDocument(self.obj.Document)
        guidoc.Document.recompute()
        guidoc.resetEdit()
        self.output_obj_mat_param()

    def output_obj_mat_param(self):
        self.print_mat_dict(self.obj.Material)
        self.print_mat_dict(self.obj.Reinforcement)
        print('\n')

    def print_mat_dict(self, mat_dict):
        if 'Name' in mat_dict:
            print('Material: {}'.format(mat_dict['Name']))
        else:
            print('Matrix material: no Name')
        for key in mat_dict:
            print('    {}: {}'.format(key, mat_dict[key]))

    # choose material card ***********************************************************************
    def get_material_card(self, material):
        for a_mat in self.materials:
            unmatched_items = set(self.materials[a_mat].items()) ^ set(material.items())
            # print(a_mat + '  -->  unmatched_items = ' + str(len(unmatched_items)))
            if len(unmatched_items) == 0:
                return a_mat
        return ""

    def choose_material_m(self, index):
        if index < 0:
            return
        # get the whole card path
        self.card_path_m = self.parameterWidget.cb_materials_m.itemData(index)
        # print('choose_material: ' + self.card_path_m)
        self.matrix = self.materials[self.card_path_m]
        self.parameterWidget.cb_materials_m.setCurrentIndex(index)
        gen_mat_desc = ""
        gen_mat_name = ""
        if 'Description' in self.matrix:
            gen_mat_desc = self.matrix['Description']
        if 'Name' in self.matrix:
            gen_mat_name = self.matrix['Name']
        self.parameterWidget.l_description_m.setText(gen_mat_desc)
        self.parameterWidget.l_name_m.setText(gen_mat_name)

    def choose_material_r(self, index):
        if index < 0:
            return
        # get the whole card path
        self.card_path_r = self.parameterWidget.cb_materials_r.itemData(index)
        # print('choose_material: ' + self.card_path_r)
        self.reinforcement = self.materials[self.card_path_r]
        self.parameterWidget.cb_materials_r.setCurrentIndex(index)
        gen_mat_desc = ""
        gen_mat_name = ""
        if 'Description' in self.reinforcement:
            gen_mat_desc = self.reinforcement['Description']
        if 'Name' in self.reinforcement:
            gen_mat_name = self.reinforcement['Name']
        self.parameterWidget.l_description_r.setText(gen_mat_desc)
        self.parameterWidget.l_name_r.setText(gen_mat_name)

    # transient material is needed if the user changed mat parameter by the mat editor
    def set_transient_material_m(self):
        self.card_path_m = '_Transient_Matrix_Material'
        self.materials[self.card_path_m] = self.matrix  # = the current matrix mat dict
        index = self.parameterWidget.cb_materials_m.findData(self.card_path_m)
        self.choose_material_m(index)

    def add_transient_material_m(self):
        self.has_transient_mat_m = True
        self.card_path_m = '_Transient_Matrix_Material'
        self.parameterWidget.cb_materials_m.addItem(
            QtGui.QIcon(":/icons/help-browser.svg"),
            self.card_path_m,
            self.card_path_m
        )
        self.set_transient_material_m()

    def set_transient_material_r(self):
        self.card_path_r = '_Transient_Reinforcement_Material'
        self.materials[self.card_path_r] = self.reinforcement  # = the current reinforced mat dict
        index = self.parameterWidget.cb_materials_r.findData(self.card_path_r)
        self.choose_material_r(index)

    def add_transient_material_r(self):
        self.has_transient_mat_r = True
        self.card_path_r = '_Transient_Reinforcement_Material'
        self.parameterWidget.cb_materials_r.addItem(
            QtGui.QIcon(":/icons/help-browser.svg"),
            self.card_path_r,
            self.card_path_r
        )
        self.set_transient_material_r()

    # edit material parameter ********************************************************************
    # TODO, also all mat parameter checks should be moved to material editor
    # and mat parameter checks should be done on analysis precheck in according to the analysis
    # should be checked if all needed parameter are defined and have all right values and units
    # TODO, if a mat card is chosen in mat editor, even without any change
    # a transient material is created too
    def edit_material_m(self):
        # opens the material editor to choose a material or edit material params
        new_material_params = self.matrix.copy()
        new_material_params = MaterialEditor.editMaterial(new_material_params)
        # material editor returns the mat_dict only, not a card_path
        # if the material editor was canceled a empty params dict will be returned,
        # do not change the self.matrix
        if new_material_params:
            self.print_mat_dict(new_material_params)
            self.matrix = new_material_params
            self.card_path_m = self.get_material_card(self.matrix)
            # TODO: mat editor returns lizence and card_name and thus the dict is not equal
            # thus the card is not found
            print('card_path: ' + self.card_path_m)
            if not self.card_path_m:
                if self.has_transient_mat_m is False:
                    self.add_transient_material_m()
                else:
                    self.set_transient_material_m()
            else:
                # we found our exact material in self.materials dict :-)
                FreeCAD.Console.PrintMessage(
                    "Material card was found in material directories. "
                    "We will use this material.\n"
                )
                index = self.parameterWidget.cb_materials_m.findData(self.card_path_m)
                # print(index)
                # set the current material in the cb widget
                self.choose_material_m(index)
        else:
            FreeCAD.Console.PrintMessage('No changes where made by the material editor.\n')

    def edit_material_r(self):
        # opens the material editor to choose a material or edit material params
        new_material_params = self.reinforcement.copy()
        new_material_params = MaterialEditor.editMaterial(new_material_params)
        # material editor returns the mat_dict only, not a card_path
        # if the material editor was canceled a empty params dict will be returned,
        # do not change the self.reinforcement
        if new_material_params:
            self.print_mat_dict(new_material_params)
            self.reinforcement = new_material_params
            self.card_path_r = self.get_material_card(self.reinforcement)
            # TODO: mat editor returns lizence and card_name and thus the dict is not equal
            # thus the card is not found
            print('card_path: ' + self.card_path_r)
            if not self.card_path_r:
                if self.has_transient_mat_r is False:
                    self.add_transient_material_r()
                else:
                    self.set_transient_material_r()
            else:
                # we found our exact material in self.materials dict :-)
                FreeCAD.Console.PrintMessage(
                    "Material card was found in material directories. "
                    "We will use this material.\n"
                )
                index = self.parameterWidget.cb_materials_r.findData(self.card_path_r)
                # print(index)
                # set the current material in the cb widget
                self.choose_material_r(index)
        else:
            FreeCAD.Console.PrintMessage('No changes where made by the material editor.\n')

    # material import ****************************************************************************
    # duplicate of std material task panel
    # TODO get rid of these redundance
    # some could be moved to material module
    def import_materials(self):
        self.pathList = []
        self.parameterWidget.cb_materials_m.clear()
        self.parameterWidget.cb_materials_r.clear()

        self.fem_prefs = FreeCAD.ParamGet(
            "User parameter:BaseApp/Preferences/Mod/Material/Resources"
        )
        self.import_solid_materials()
        # self.print_materialsdict()

    def import_solid_materials(self):
        use_built_in_materials = self.fem_prefs.GetBool("UseBuiltInMaterials", True)
        if use_built_in_materials:
            system_mat_dir = FreeCAD.getResourceDir() + "/Mod/Material/StandardMaterial"
            self.add_cards_from_a_dir(system_mat_dir, ":/icons/freecad.svg")

        use_mat_from_config_dir = self.fem_prefs.GetBool("UseMaterialsFromConfigDir", True)
        if use_mat_from_config_dir:
            user_mat_dirname = FreeCAD.getUserAppDataDir() + "Material"
            self.add_cards_from_a_dir(user_mat_dirname, ":/icons/preferences-general.svg")

        use_mat_from_custom_dir = self.fem_prefs.GetBool("UseMaterialsFromCustomDir", True)
        if use_mat_from_custom_dir:
            custom_mat_dir = self.fem_prefs.GetString("CustomMaterialsDir", "")
            self.add_cards_from_a_dir(custom_mat_dir, ":/icons/user.svg")

    def add_cards_from_a_dir(self, mat_dir, icon):
        # fill list of card names and
        import glob
        from importFCMat import read
        dir_path_list = glob.glob(mat_dir + '/*' + ".FCMat")
        self.pathList = self.pathList + dir_path_list

        # fill self.materials
        for a_path in dir_path_list:
            mat_dict = read(a_path)
            # check if the dict exists in materials
            # TODO if the unit is different two cards would be different too
            if mat_dict not in self.materials.values():
                self.materials[a_path] = mat_dict
                self.icons[a_path] = icon

    def add_cards_to_combo_boxes(self):
        # fill comboboxes, in combo boxes the card name is used not the material name
        import os
        card_name_list = []
        for a_path in self.materials:
            card_name = os.path.basename(a_path[:-(len(".FCMat"))])
            card_name_list.append([card_name, a_path, self.icons[a_path]])
        card_name_list.sort()
        for mat in card_name_list:
            self.parameterWidget.cb_materials_m.addItem(QtGui.QIcon(mat[2]), mat[0], mat[1])
        for mat in card_name_list:
            self.parameterWidget.cb_materials_r.addItem(QtGui.QIcon(mat[2]), mat[0], mat[1])
