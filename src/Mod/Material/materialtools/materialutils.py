# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2018 - Bernd Hahnebach <bernd@bimstatik.org>            *
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



def getMaterialResources():
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
        resources.append(FreeCAD.getResourceDir() + "Mod" + os.sep + "Material" + os.sep + "StandardMaterial")
    if use_mat_from_config_dir:
        resources.append(FreeCAD.ConfigGet("UserAppData") + "Material")
    if use_mat_from_custom_dir:
        custom_mat_dir = fem_prefs.GetString("CustomMaterialsDir", "")
        if os.path.exists(custom_mat_dir):
            resources.append(custom_mat_dir)
    outputResources(resources)
    return resources


def outputResources(resources):
    print('locations we gone look for material cards:')
    for path in resources:
        print('  ' + path)
    print('\n')


def outputCards(cards):
    print('material cards:')
    for card in cards:
        print('  ' + card + ': ' + cards[card])
    print('\n')


def getMaterialCards(resources):
    # no duplicates, an existing card will be overwritten!
    cards = {}
    for p in resources:
        for f in os.listdir(p):
            b, e = os.path.splitext(f)
            if e.upper() == ".FCMAT":
                # print(type(b))
                if isinstance(b, str):  # checks for string, returns false for a unicode string
                    b = b.decode('utf-8')  # qt needs unicode to display the special characters the right way
                cards[b] = p + os.sep + f
    # outputCards()
    getAndOutputAllCardData()
    # print(cards)
    return cards


def getAndOutputAllCardData(cards):
    print('\n\n\nMYSTART')
    # get all registered material property keys
    registedCardKeys = []
    for gr in getMaterialAttributeStructure():
        for key in gr[1]:
            registedCardKeys.append(key)
    registedCardKeys = sorted(registedCardKeys)

    # get all data from all known cards
    allCardsAndData = {}  # {cardfilename: ['path', materialdict]}
    for card in cards:
        import importFCMat
        d = importFCMat.read(cards[card])
        allCardsAndData[card] = [cards[card], d]
    '''
    for card in allCardsAndData:
        print(card)
        print(allCardsAndData[card][0])
        print(allCardsAndData[card][1])
        print('\n')
    '''

    # find not registered and registered keys in the used data
    usedAndRegisteredCardKeys = []
    usedAndNotRegisteredCardKeys = []
    registeredAndNotUsedCardKeys = []
    for card in allCardsAndData:
        for k in allCardsAndData[card][1]:
            if k in registedCardKeys:
                usedAndRegisteredCardKeys.append(k)
            else:
                usedAndNotRegisteredCardKeys.append(k)
    for k in registedCardKeys:
        if (k not in usedAndRegisteredCardKeys) and (k not in usedAndNotRegisteredCardKeys):
            registeredAndNotUsedCardKeys.append(k)
    usedAndRegisteredCardKeys = sorted(list(set(usedAndRegisteredCardKeys)))
    usedAndNotRegisteredCardKeys = sorted(list(set(usedAndNotRegisteredCardKeys)))
    registeredAndNotUsedCardKeys = sorted(list(set(registeredAndNotUsedCardKeys)))
    print(usedAndRegisteredCardKeys)
    print(usedAndNotRegisteredCardKeys)
    print(registeredAndNotUsedCardKeys)
    # still there may be lots of properties in the template which are not used in other materials but the tmplate is handeled here like a material
    print('MYEND\n\n\n')
