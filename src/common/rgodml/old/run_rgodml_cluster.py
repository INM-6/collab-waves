# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:52:12 2014

@author: zehl
"""
import os

import ox.factory as oxf

from rgodml import RG_sections
from rgodml import RG_properties
from rgodml import RG_values
from rgodml import blackrock_templates

from rgcollector import RG_FileCollector

author = "Lyuba Zehl"

monkey = "Lilou"
fileformat = "odml"

odml_directory = os.path.join(os.sep, 'scratch', 'zehl', 'MetaData'+monkey, 'odMLbuild', '')
data_directory = os.path.join(os.sep, 'datasets', 'marseille', 'congloue', 'data', 'DataGrasp', '')
metadata_directory = os.path.join(os.sep, 'scratch', 'zehl', '')

# get PBS_ARRAYID of
num = int(os.getenv('PBS_ARRAYID'))

subsessions = RG_FileCollector.get_all_subsessions(monkey, data_directory+'Data'+monkey+os.sep)
subsession = subsessions[num-1]
print "Subsession:", subsession

sec_dict = RG_sections.section_dictionary()
prop_dict = RG_properties.property_dictionary()
val_data_dict = RG_values.value_data_dictionary(monkey, subsession, data_directory, metadata_directory)

blr_sec_dict, blr_prop_dict = blackrock_templates.blr_dictionary_template()

odmlobj = oxf.Factory(author = author)
odmlobj.add_section_dict(sec_dict)
odmlobj.add_property_dict(prop_dict)

odmlobj.add_section_dict(blr_sec_dict)
odmlobj.add_property_dict(blr_prop_dict)

for prop_path in val_data_dict.keys():
    odmlobj.replace_value_data(val_data_dict[prop_path], prop_path)

odmlobj.save(odml_directory, subsession, fileformat = fileformat)