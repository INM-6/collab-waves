# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 22:49:29 2013

@author: zehl
"""
import os
import sys
import traceback

import ox.factory as oxf

import rgodml.RG_Sections as rgs
import rgodml.RG_Properties as rgp
import rgodml.RG_Values as rgv

from rgcollector import RG_FileCollector as fcol

#==============================================================================
# Define parameters needed for running the odml generation.
# monkey: given name of the monkey ('Lilou', 'Tanya', 'Tanya2' or 'Nikos')
# subsessions: list of subsession names (e.g.: ['l101015-001', 'l101013-002'])
# datadir: directory to the recorded data
# mdatadir: directory to the metadata
# odmldir: directory to the odml
#==============================================================================
if len(sys.argv) < 5 or len(sys.argv) > 6:
    raise ValueError(
          "Number of arguments in the command line option is wrong. \n "+
          "**Args:** \n script_location  \n monkey \n subsession \n "+
          "author (firstname-lastname) \n startfrom_subsession")

if sys.argv[2] in ['Lilou', 'Nikos', 'Tanya', 'Tanya2']:
    monkey = sys.argv[2]

else:
    raise ValueError(sys.argv[2] +
          " is not valid as 2nd argument in the command line option. \n "+
          "Valid Args: 'Lilou', 'Nikos', 'Tanya' or 'Tanya2'")


if sys.argv[1] == 'marseille':
    main_dir = os.path.join(os.sep, 'riou', 'work', 'comco_congloue',
                            'DataGrasp')
    data_dir = os.path.join(main_dir)

    mdata_dir = os.path.join(main_dir)

    odml_dir = os.path.join(main_dir, 'MetaData'+monkey, 'odMLbuild')

elif sys.argv[1] == 'juelich':
    main_dir = os.path.join(os.sep, 'media', 'Backup', 'zehl', 'datasets',
                            'DataGrasp')
    data_dir = os.path.join(main_dir)

    mdata_dir = os.path.join(os.sep, 'home', 'zehl', 'datasets', 'marseille',
                             'DataGrasp')
    odml_dir = os.path.join(os.sep, 'home', 'zehl', 'datasets', 'marseille',
                             'DataGrasp', 'MetaData'+monkey, 'odMLbuild')

else:
    raise ValueError(sys.argv[1] +
          " is not a valid argument for the command line option. \n "+
          "Valid Args: 'marseille' or 'juelich'")

if sys.argv[3] == 'all':
    subsessions = fcol.get_all_subsessions(monkey, main_dir)

else:
    subsessions = [sys.argv[3]]


author = ' '.join(sys.argv[4].split('-'))

#==============================================================================
# Redifine subsessions list if you want to start from a certain subsession
#==============================================================================
try:
    startfrom_subsession = sys.argv[5]
    subsessions = subsessions[subsessions.index(startfrom_subsession):]
except:
    pass

#==============================================================================
# run odml generation on LAPTOP
#=============================================================================+

for num in range(len(subsessions)):
    subsession = subsessions[num]
#    print subsession, num+1, len(subsessions)
    try:

        sec_dict = rgs.section_dictionary()
        prop_dict = rgp.property_dictionary()
        val_dict1 = rgv.value_dictionary(monkey, subsession, data_dir, mdata_dir)
        sort_ver, val_dict2 = rgv.spikesorting_value_dictionary(monkey, subsession, data_dir, mdata_dir)

        fobj = oxf.Factory(author)
        fobj.add_section_dict(sec_dict)
        fobj.add_property_dict(prop_dict)
        fobj.add_value_dict(val_dict1)

        if sort_ver:
            fobj.duplicate_section('/SpikeSorting/VersionX', sort_ver)
            fobj.add_value_dict(val_dict2)
        else: pass

        fobj.save(fobj.temp, odml_dir, subsession, 'odml')

    except Exception:
        print("Subsession %s does not work. The following Error occurred:") %(subsession)
        print(traceback.print_exc())
        print
        pass