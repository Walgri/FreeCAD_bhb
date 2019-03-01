# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2013-2015 - Juergen Riegel <FreeCAD@juergen-riegel.net> *
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


import FreeCAD


# here the usage description if you use this tool from the command line ("__main__")
CommandlineUsage = """Material - Tool to work with FreeCAD Material definition cards

Usage:
   Material [Options] card-file-name

Options:
 -c, --output-csv=file-name     write a comma separated grid with the material data

Exit:
 0      No Error or Warning found
 1      Argument error, wrong or less Arguments given

Tool to work with FreeCAD Material definition cards

Examples:

   Material  "StandardMaterial/Steel.FCMat"

Author:
  (c) 2013 Juergen Riegel
  mail@juergen-riegel.net
  Licence: LGPL

Version:
  0.1
"""


# see comments in module importFCMat, there is an independent parser implementation for reading and writing FCMat files
# inside FreeCAD a mixture of these parsers and the ones in importFCMat.py is used


def importFCMat(fileName):
    "Read a FCMat file into a dictionary"
    try:
        import ConfigParser as configparser
    except ImportError:
        import configparser

    Config = configparser.RawConfigParser()
    Config.optionxform = str
    Config.read(fileName)
    dict1 = {}
    for section in Config.sections():
        options = Config.options(section)
        for option in options:
            dict1[option] = Config.get(section, option)

    return dict1


def exportFCMat(fileName, matDict):
    "Write a material dictionary to a FCMat file"
    try:
        import ConfigParser as configparser
    except ImportError:
        import configparser
    import string
    Config = configparser.RawConfigParser()

    # create groups
    for x in matDict.keys():
        grp, key = string.split(x, sep='_')
        if not Config.has_section(grp):
            Config.add_section(grp)

    # fill groups
    for x in matDict.keys():
        grp, key = string.split(x, sep='_')
        Config.set(grp, key, matDict[x])

    Preamble = "# This is a FreeCAD material-card file\n\n"
    # Writing our configuration file to 'example.cfg'
    with open(fileName, 'wb') as configfile:
        configfile.write(Preamble)
        Config.write(configfile)


def get_material_template(withSpaces=False):
    # material properties
    # see the following resources in the FreeCAD wiki for more information about the material specific properties:
    # https://www.freecadweb.org/wiki/Material_data_model
    # https://www.freecadweb.org/wiki/Material

    import yaml
    template_data = yaml.safe_load(open(FreeCAD.ConfigGet('AppHomePath') + 'Mod/Material/Templatematerial.yml'))
    if withSpaces:
        # on attributes, add a space before a capital letter, will be used for better display in the ui
        import re
        for group in template_data:
            gg = list(group.keys())[0]  # group dict has only one key
            for proper in list(group[gg].keys()):  # iterating over a dict and changing it is not allowed, thus we iterate over a list of the keys
                new_proper = re.sub(r"(\w)([A-Z]+)", r"\1 \2", proper)
                group[gg][new_proper] = group[gg][proper]
                del group[gg][proper]
    return template_data


def create_mat_tools_header():
    template_data = get_material_template()
    # tools will not be copied ... thus they are not there ...
    # TODO copy them if they are reimplemented in Python
    # headers = FreeCAD.ConfigGet('AppHomePath') + 'Mod/Material/StandardMaterial/Tools/headersnew'
    headers = '/home/hugo/Documents/dev/freecad/freecadbhb_dev/freecad/src/Mod/Material/StandardMaterial/Tools/headers'
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
    rev = FreeCAD.ConfigGet("BuildVersionMajor") + "." + FreeCAD.ConfigGet("BuildVersionMinor") + "." + FreeCAD.ConfigGet("BuildRevision")
    template_data = get_material_template()
    template_card = '/home/hugo/Documents/dev/freecad/freecadbhb_dev/freecad/src/Mod/Material/StandardMaterial/TEMPLATE.FCMat'
    f = open(template_card, "w")
    f.write('; TEMPLATE\n')
    f.write('; (c) 2013-2015 Juergen Riegel (CC-BY 3.0)\n')
    f.write('; information about the content of such cards can be found on the wiki:\n')
    f.write('; https://www.freecadweb.org/wiki/index.php?title=Material\n')
    f.write(': this template card was created by FreeCAD ' + rev + '\n\n')
    f.write('; localized Name, Description and KindOfMaterial uses 2 letter codes\n')
    f.write('; defined in ISO-639-1, see https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes\n')
    f.write('; find unit information in src/App/FreeCADInit.py\n')
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
                f.write("\n[" + gg + "]\n")
            for prop_name in group[gg]:
                description = group[gg][prop_name]['Description']
                if not description.strip():
                    f.write('; Description to be updated\n')
                else:
                    f.write('; ' + description + '\n')
                url = group[gg][prop_name]['URL']
                if url.strip():
                    f.write('; ' + url + '\n')
                f.write(prop_name + ' =\n\n')
    f.close


def read_cards_from_path(cards_path):
    from os import listdir
    from os.path import isfile, join, basename, splitext
    from importFCMat import read
    only_files = [f for f in listdir(cards_path) if isfile(join(cards_path, f))]
    mat_files = [f for f in only_files if basename(splitext(f)[1]) == '.FCMat' or basename(splitext(f)[1]) == '.fcmat']
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


if __name__ == '__main__':
    import sys
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:", ["output-csv="])
    except getopt.GetoptError:
        # print help information and exit:
        sys.stderr.write(CommandlineUsage)
        sys.exit(1)

    # checking on the options
    for o, a in opts:
        if o in ("-c", "--output-csv"):
            print("writing file: " + a + "\n")
            OutPath = a

    # running through the files
    FileName = args[0]

    kv_map = importFCMat(FileName)
    for k in kv_map.keys():
        print(repr(k) + " : " + repr(kv_map[k]))
    sys.exit(0)  # no error
