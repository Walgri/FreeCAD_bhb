# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2019 - Bernd Hahnebach <bernd@bimstatik.org>            *
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

__title__ = "material utilities"
__author__ = "Bernd Hahnebach"
__url__ = "http://www.freecadweb.org"

import os
import sys
from os.path import join

import FreeCAD

if sys.version_info.major >= 3:
    unicode = str


# card tools will not be copied inot build ... thus they are not there ...
# TODO copy them, than setting users src path is no longer needed
# headers = FreeCAD.ConfigGet('AppHomePath') + 'Mod/Material/StandardMaterial/Tools/headers'
FREECAD_SRC_PATH = '/home/hugo/Documents/dev/freecad/freecadbhb_dev/freecad/src'


# TODO: may be move the module into a package materialtools/utils, do not forget the __init__.py
# cmake is much easier without a package, thus leave it for now in main module directory
# TODO: rename to cardutils or cardtools if inside materialtools


def get_material_resources():
    fem_prefs = FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Mod/Material/Resources")
    use_built_in_materials = fem_prefs.GetBool("UseBuiltInMaterials", True)
    use_mat_from_config_dir = fem_prefs.GetBool("UseMaterialsFromConfigDir", True)
    use_mat_from_custom_dir = fem_prefs.GetBool("UseMaterialsFromCustomDir", True)
    if use_mat_from_custom_dir:
        custom_mat_dir = fem_prefs.GetString("CustomMaterialsDir", "")
    # later found cards with same name will override cards
    # FreeCAD returns paths with / at the end, thus not os.sep is needed on first +
    resources = []
    if use_built_in_materials:
        resources.append(
            join(FreeCAD.getResourceDir(), "Mod", "Material", "StandardMaterial")
        )
    if use_mat_from_config_dir:
        resources.append(FreeCAD.ConfigGet("UserAppData") + "Material")
    if use_mat_from_custom_dir:
        custom_mat_dir = fem_prefs.GetString("CustomMaterialsDir", "")
        if os.path.exists(custom_mat_dir):
            resources.append(custom_mat_dir)
    return resources


def output_resources(resources):
    print('locations we gone look for material cards:')
    for path in resources:
        print('  ' + path)
    print('\n')


def output_cards(cards):
    print('material cards:')
    for card in cards:
        print('  ' + card + ': ' + cards[card])
    print('\n')


def get_material_cards(resources=None):

    # get resources if not given
    if not resources:
        resources = get_material_resources()
    # output_resources(resources)

    cards = {}
    # no duplicates in card names, an existing card will be overwritten!
    # means the last found is taken, if two cards with the same name exists
    for p in resources:
        if os.path.exists(p):
            for f in sorted(os.listdir(p)):  # sort cards, TODO: use pref
                b, e = os.path.splitext(f)
                if e.upper() == ".FCMAT":
                    # TODO: get special characters working in any case
                    '''
                    # print(type(b))
                    # checks for a string, returns false for a unicode string
                    if isinstance(b, str):
                        # qt needs unicode to display the special characters the right way
                        b = b.decode('utf-8')
                    '''
                    cards[b] = p + os.sep + f
    # outputCards()
    # get_and_output_all_carddata(cards)  # to compare with cards, has some problems ...
    return cards


def get_and_output_all_carddata(cards):
    print('\n\n\nMYSTART\n----------')
    # get all registered material property keys
    registed_cardkeys = []
    template_data = get_material_template()
    # print(template_data)
    for group in template_data:
        gg = list(group.keys())[0]  # group dict has only one key
        for key in group[gg]:
            registed_cardkeys.append(key)
    registed_cardkeys = sorted(registed_cardkeys)
    # print(registed_cardkeys)

    # get all data from all known cards
    all_cards_and_data = {}  # {cardfilename: ['path', materialdict]}
    for card in cards:
        from importFCMat import read
        d = read(cards[card])
        all_cards_and_data[card] = [cards[card], d]
    '''
    for card in all_cards_and_data:
        print(card)
        print(all_cards_and_data[card][0])
        print(all_cards_and_data[card][1])
        print('\n')
    '''

    # find not registered and registered keys in the used data
    used_and_registered_cardkeys = []
    used_and_not_registered_cardkeys = []
    registered_and_not_used_cardkeys = []
    for card in all_cards_and_data:
        for k in all_cards_and_data[card][1]:
            if k in registed_cardkeys:
                used_and_registered_cardkeys.append(k)
            else:
                used_and_not_registered_cardkeys.append(k)
    for k in registed_cardkeys:
        if (k not in used_and_registered_cardkeys) and (k not in used_and_not_registered_cardkeys):
            registered_and_not_used_cardkeys.append(k)

    used_and_registered_cardkeys = sorted(list(set(used_and_registered_cardkeys)))
    used_and_not_registered_cardkeys = sorted(list(set(used_and_not_registered_cardkeys)))
    registered_and_not_used_cardkeys = sorted(list(set(registered_and_not_used_cardkeys)))
    FreeCAD.Console.PrintMessage(
        '\nused_and_registered_cardkeys:\n{}\n'
        .format(used_and_registered_cardkeys)
    )
    FreeCAD.Console.PrintMessage(
        '\nused_and_not_registered_cardkeys:\n{}\n'
        .format(used_and_not_registered_cardkeys)
    )
    FreeCAD.Console.PrintMessage(
        '\nregistered_and_not_used_cardkeys:\n{}\n'
        .format(registered_and_not_used_cardkeys)
    )

    # still there might be lots of properties in the template
    # which are not used in other materials
    # but the tmplate is handeled here like a material
    print('----------\nMYEND\n\n\n')


def create_mat_tools_header():
    template_data = get_material_template()
    headers = join(FREECAD_SRC_PATH, 'Mod/Material/StandardMaterial/Tools/headers')
    f = open(headers, "w")
    for group in template_data:
        gg = list(group.keys())[0]  # group dict has only one key
        # do not write group UserDefined
        if gg != 'UserDefined':
            for prop_name in group[gg]:
                if prop_name != 'None':
                    f.write(prop_name + '\n')
    f.close


def create_mat_template_card(write_group_section=True):
    rev = "{}.{}.{}".format(
        FreeCAD.ConfigGet("BuildVersionMajor"),
        FreeCAD.ConfigGet("BuildVersionMinor"),
        FreeCAD.ConfigGet("BuildRevision")
    )
    template_data = get_material_template()
    template_card = join(FREECAD_SRC_PATH, 'Mod/Material/StandardMaterial/TEMPLATE.FCMat')
    f = open(template_card, "w")
    f.write('; TEMPLATE\n')
    f.write('; (c) 2013-2015 Juergen Riegel (CC-BY 3.0)\n')
    f.write('; information about the content of such cards can be found on the wiki:\n')
    f.write('; https://www.freecadweb.org/wiki/Material\n')
    f.write(': this template card was created by FreeCAD ' + rev + '\n\n')
    f.write('; localized Name, Description and KindOfMaterial uses 2 letter codes\n')
    f.write('; defined in ISO-639-1, see https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes\n')
    f.write('; find unit information in src/App/FreeCADInit.py')
    # write sections
    # write standard FCMat section if write group section parameter is set to False
    if write_group_section is False:
        f.write("\n[FCMat]\n")
    for group in template_data:
        gg = list(group.keys())[0]  # group dict has only one key
        # do not write groups Meta and UserDefined
        if (gg != 'Meta') and (gg != 'UserDefined'):
            # only write group section if write group section parameter is set to True
            if write_group_section is True:
                f.write("\n\n[" + gg + "]")
            for prop_name in group[gg]:
                f.write('\n')
                description = group[gg][prop_name]['Description']
                if not description.strip():
                    f.write('; Description to be updated\n')
                else:
                    f.write('; ' + description + '\n')
                url = group[gg][prop_name]['URL']
                if url.strip():
                    f.write('; ' + url + '\n')
                f.write(prop_name + ' =\n')
    f.close


def get_material_template(withSpaces=False):
    # material properties
    # see the following resources in the FreeCAD wiki for more
    # information about the material specific properties:
    # https://www.freecadweb.org/wiki/Material_data_model
    # https://www.freecadweb.org/wiki/Material

    import yaml
    template_data = yaml.safe_load(
        open(join(FreeCAD.ConfigGet('AppHomePath'), 'Mod/Material/Templatematerial.yml'))
    )
    if withSpaces:
        # on attributes, add a space before a capital letter
        # will be used for better display in the ui
        import re
        for group in template_data:
            gg = list(group.keys())[0]  # group dict has only one key
            # iterating over a dict and changing it is not allowed
            # thus it is iterated over a list of the keys
            for proper in list(group[gg].keys()):
                new_proper = re.sub(r"(\w)([A-Z]+)", r"\1 \2", proper)
                group[gg][new_proper] = group[gg][proper]
                del group[gg][proper]
    return template_data


def read_cards_from_path(cards_path):
    from os import listdir
    from os.path import isfile, join, basename, splitext
    from importFCMat import read
    only_files = [f for f in listdir(cards_path) if isfile(join(cards_path, f))]
    mat_files = [f for f in only_files if basename(splitext(f)[1]) == '.FCMat']
    # print(mat_files)
    mat_cards = []
    for f in sorted(mat_files):
        mat_cards.append(read(join(cards_path, f)))
    return mat_cards


def write_cards_to_path(cards_path, cards_data, write_group_section=True, write_template=False):
    from importFCMat import write
    from os.path import join
    for card_data in cards_data:
        if (card_data['CardName'] == 'TEMPLATE') and (write_template is False):
            continue
        else:
            card_path = join(cards_path, (card_data['CardName'] + '.FCMat'))
            print(card_path)
            if write_group_section is True:
                write(card_path, card_data, True)
            else:
                write(card_path, card_data, False)
