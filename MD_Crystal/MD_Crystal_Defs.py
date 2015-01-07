# coding=UTF-8
"""
  MD_Crystal_Defs.py
  
  Authors:  Ardiani, Franco
            Manelli, Andrés
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
  MA 02110-1301, USA.
  
  ------------------------------------------------------------------

  Python Module MD_Crystal
  
"""

from ovito.modifiers import BondAngleAnalysisModifier

OTHER = BondAngleAnalysisModifier.Type.OTHER
FCC = BondAngleAnalysisModifier.Type.FCC
HCP = BondAngleAnalysisModifier.Type.HCP
BCC = BondAngleAnalysisModifier.Type.BCC
ICO = BondAngleAnalysisModifier.Type.ICO

STRUCT_TYPES = [OTHER, FCC, HCP, BCC, ICO]
