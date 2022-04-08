"""
Script to add immunogenicity annotations to NeoDisc files. The immunogenicity annotations are extracted from the
supplementary files of the respective publications. The names and location of these files are configured in the
Parameters class.

@Author: Markus Muller, CHUV/Ludwig, Lausanne, Switzerland
@License : (C)Copyright 2022
"""

from DataWrangling.RosenbergImmunogenicityAnnotatorShort import *
from DataWrangling.RosenbergImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorLong import *
from DataWrangling.NeoDiscImmunogenicityAnnotatorShort import *
from DataWrangling.TESLAImmunogenicityAnnotatorLong import *
from DataWrangling.TESLAImmunogenicityAnnotatorShort import *
from Utils.DataManager import *

# load valid patients (don't check for presence of immunogenic peptides yet)
dataManager = DataManager(immunogenity=False)

converter = RosenbergImmunogenicityAnnotatorLong(dataManager)
converter.annotate_response_types()

converter = RosenbergImmunogenicityAnnotatorShort(dataManager)
converter.annotate_response_types()

#converter = NeoDiscImmunogenicityAnnotatorShort(dataManager)
#converter.annotate_response_types()

# converter = NeoDiscImmunogenicityAnnotatorLong(dataManager)
# converter.annotate_response_types()

#converter = TESLAImmunogenicityAnnotatorLong(dataManager)
#converter.annotate_response_types()

#converter = TESLAImmunogenicityAnnotatorShort(dataManager)
#converter.annotate_response_types()

dataManager = DataManager(immunogenity=True)
