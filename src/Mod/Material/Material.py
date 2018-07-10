#***************************************************************************
#*                                                                         *
#*   Copyright (c) 2013-2015 - Juergen Riegel <FreeCAD@juergen-riegel.net> *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU Lesser General Public License (LGPL)    *
#*   as published by the Free Software Foundation; either version 2 of     *
#*   the License, or (at your option) any later version.                   *
#*   for detail see the LICENCE text file.                                 *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU Library General Public License for more details.                  *
#*                                                                         *
#*   You should have received a copy of the GNU Library General Public     *
#*   License along with this program; if not, write to the Free Software   *
#*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
#*   USA                                                                   *
#*                                                                         *
#***************************************************************************


def getMaterialAttributeStructure(withSpaces=None):
    # material properties
    # see the following resources in the FreeCAD wiki for more information about the material specific properties:
    # https://www.freecadweb.org/wiki/Material_data_model
    # https://www.freecadweb.org/wiki/Material
    materialPropertyGroups = (
        ("Meta", (
            "CardName",
            "AuthorAndLicense",
            "Source"
        )),
        ("General", (
            "Name",
            "Father",
            "Description",
            "Density",
            "Vendor",
            "ProductURL",
            "SpecificPrice"
        )),
        ("Mechanical", (
            "YoungsModulus",  # https://en.wikipedia.org/wiki/Young%27s_modulus
            "PoissonRatio",  # https://en.wikipedia.org/wiki/Poisson%27s_ratio
            "UltimateTensileStrength",  # https://en.wikipedia.org/wiki/Ultimate_tensile_strength
            "CompressiveStrength",  # https://en.wikipedia.org/wiki/Compressive_strength
            "YieldStrength",  # https://en.wikipedia.org/wiki/Yield_Strength
            "UltimateStrain",  # https://en.wikipedia.org/wiki/Ultimate_tensile_strength
            "FractureToughness",  # https://en.wikipedia.org/wiki/Fracture_toughness
            "AngleOfFriction"  # https://en.wikipedia.org/wiki/Friction#Angle_of_friction and https://en.m.wikipedia.org/wiki/Mohr%E2%80%93Coulomb_theory
        )),
        ("Thermal", (
            "ThermalConductivity",  # https://en.wikipedia.org/wiki/Thermal_conductivity
            "ThermalExpansionCoefficient",  # https://en.wikipedia.org/wiki/Volumetric_thermal_expansion_coefficient
            "SpecificHeat"  # https://en.wikipedia.org/wiki/Heat_capacity
        )),
        ("Architectural", (
            "Model",
            "ExecutionInstructions",
            "FireResistanceClass",
            "StandardCode",
            "SoundTransmissionClass",
            "Color",
            "Finish",
            "UnitsPerQuantity",
            "EnvironmentalEfficiencyClass"
        )),
        ("Rendering", (
            "DiffuseColor",
            "AmbientColor",
            "SpecularColor",
            "Shininess",
            "EmissiveColor",
            "Transparency",
            "VertexShader",
            "FragmentShader",
            "TexturePath",
            "TextureScaling"
        )),
        ("Vector rendering", (
            "ViewColor",
            "ViewFillPattern",
            "SectionFillPattern",
            "ViewLinewidth",
            "SectionLinewidth"
        )),
        ("User defined", (
        ))
    )
    if withSpaces:
        # on attributes, add a space before a capital letter, will be used for better display in the ui
        import re
        newMatProp = []
        for group in materialPropertyGroups:
            newAttr = []
            for attr in group[1]:
                newAttr.append(re.sub(r"(\w)([A-Z])", r"\1 \2", attr))
            newMatProp.append([group[0], newAttr])
        materialPropertyGroups = newMatProp
    # print(materialPropertyGroups)
    return materialPropertyGroups
