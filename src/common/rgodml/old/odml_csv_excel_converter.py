# -*- coding: utf-8 -*-
"""
Created on Tue May 27 10:05:50 2014

@author: zehl
"""
import csv
import odml
import xlrd

def odmltree2csv4excel(document, save_to):
    """
    Writing an empty odml tree structure to an csv usable in Excel.
    """

    tree = [{'Section Name': s.name,
             'Section Type': s.type,
             'Path 2 Section': s.get_path(),
             'Section Definition': s.definition} \
             for s in document.itersections()]

    header = ['Path 2 Section',
              'Section Name',
              'Section Type',
              'Section Definition']

    csvfile = open(save_to, 'wb')
    csvwriter = csv.DictWriter(csvfile, fieldnames=header, dialect='excel')

    csvwriter.writeheader()
    csvwriter.writerows(tree)


def odmlprop2csv4excel(document, save_to):
    """
    Writing odml properties of an odml tree to an csv usable in Excel.
    """

    tree = [{'Property Name': v.parent.name,
             'Property Definition': v.parent.definition,
             'Path 2 Property': v.parent.get_path(),
             'Value Definition': v.definition,
             'Value': v.data,
             'odML Data Type': v.dtype,
             'Data Unit': v.unit,
             'Data Uncertainty': v.uncertainty} \
             for v in document.itervalues()]

    header = ['Path 2 Property',
              'Property Name',
              'Property Definition',
              'Value',
              'Value Definition',
              'odML Data Type',
              'Data Unit',
              'Data Uncertainty']

    csvfile = open(save_to, 'wb')
    csvwriter = csv.DictWriter(csvfile, fieldnames=header, dialect='excel',
                               quoting=csv.QUOTE_NONNUMERIC)

    csvwriter.writeheader()
    for i in range(len(tree)):
        if i > 0 and tree[i]['Path 2 Property'] == tree[i-1]['Path 2 Property']:
            tree[i]['Path 2 Property'] = ''
            tree[i]['Property Name'] = ''
            tree[i]['Property Definition'] = ''
        csvwriter.writerow(tree[i])


def propcsv2dict(load_from):
    """
    Writing odml properties within an csv to a list of dictionaries.
    """
    header = ['Path 2 Property',
              'Property Name',
              'Property Definition',
              'Value',
              'Value Definition',
              'odML Data Type',
              'Data Unit',
              'Data Uncertainty']

    csvfile = open(load_from)
    csvreader = csv.DictReader(csvfile, fieldnames=header, dialect='excel',
                       quoting=csv.QUOTE_NONNUMERIC)

    pdict = {}
    for i, row in enumerate(csvreader):
        if i > 0 and pdict.has_key(row['Property Name']):
            pdict[row['Property Name']]['Value'].append(row['Value'])
            pdict[row['Property Name']]['Value Definition'].append(row['Value Definition'])
            pdict[row['Property Name']]['odML Data Type'].append(row['odML Data Type'])
            pdict[row['Property Name']]['Data Unit'].append(row['Data Unit'])
            pdict[row['Property Name']]['Data Uncertainty'].append(row['Data Uncertainty'])
        elif i > 0 and pdict.has_key(row['Property Name']) == False:
            pdict[row['Property Name']] = {'Path 2 Property': row['Path 2 Property'],
                                           'Property Definition': row['Property Definition'],
                                           'Value': [row['Value']],
                                           'Value Definition': [row['Value Definition']],
                                           'odML Data Type': [row['odML Data Type']],
                                           'Data Unit': [row['Data Unit']],
                                           'Data Uncertainty': [row['Data Uncertainty']]}

    return pdict


def propcsv2odmltree(load_from, document):
    """
    Writing odml properties within an csv to an odml tree.
    """
    pdict = propcsv2dict(load_from)

    for pname in pdict.keys():
        values = []
        for i in range(len(pdict[pname]['Value'])):
            values.append(odml.Value(data=pdict[pname]['Value'][i],
                                     dtype=pdict[pname]['odML Data Type'][i],
                                     unit=pdict[pname]['Data Unit'],
                                     uncertainty=pdict[pname]['Data Uncertainty'],
                                     definition=pdict[pname]['Value Definition']))

        prop = document.get_property_by_path(pdict[pname]['Path 2 Property'])
        prop.definition = pdict[pname]['Property Definition']
        prop.values = values


def xls_sheet_reader(load_from, sheetname, num_parent_key_cols, parent_key_sep, num_header_rows):
    """
    Reads an xls (Microsoft Excel) sheet containing metadata information.

    This function reads the information on metadata contained in an Excel
    file and returns a nested dictionary containing the information. The top
    level dictionary contains one (parent) key per row which is named according
    to the concatenation first num_parent_key_cols columns in that row of the
    sheet, separated by parent_key_sep. The corresponding parent values are
    themselves nested dictionaries that replicate the structure given by the top
    num_header_rows rows of the sheet for columns greater or equal to
    num_parent_key_cols (each row is one nesting level).

    Args:
        filename (str):
            Filename of the Excel file
        sheetname (str):
            Name of the Excel sheet
        num_parent_key_cols (int):
            Number of columns from which the parent key of each row is generated
        parent_key_sep (str):
            String that separates the num_parent_key_cols components that make
            up the parent key if num_parent_key_cols>1.
        num_header_rows (int):
            Number of rows that contain the nested header structure in the Excel
            file. The top row (row 0) is the top level heading, ... below the
            parent key

    Return:
        dict:
            A nested dictionary containing the Excel file contents
        list:
            A list containing the raw contents of the Excel file, with the
            num_header_rows headers duplicated where necessary. This list may be
            passed directly to the csvwriter.writerows() method of the csv
            package to easily create a matching and compatible csv output:
                import csv
                dictout, csvout = xls_sheet_reader(filename+'.xls', 'Sheet1', 2, '', 2)
                csvbook = csv.writer(open(filename + '.csv', 'wb'), delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                csvbook.writerows(csvout)
    """

    # open excel sheet for reading
    xlsbook = xlrd.open_workbook(load_from)
    sheet = xlsbook.sheet_by_name(sheetname)

    # result structure
    result = {}

    # contains each header row
    header = []

    # header rows
    for headers in xrange(num_header_rows):
        header.append(sheet.row(headers))

        # fill empty headings on top row
        # (spanning columns)
        for heading in xrange(1, len(header[headers])):
            # empty?
            if header[headers][heading].value == u'':
                # only fill if this cell and the one to the left
                # have the same parent (or if this is the top heading)
                if headers == 0 or (headers > 0 and header[headers - 1][heading] == header[headers - 1][heading - 1]):
                    header[headers][heading] = header[headers][heading - 1]

    # add headers to csv output
    csvoutput = []
    for headers in xrange(num_header_rows):
        csvrow = []
        for heading in xrange(0, len(header[headers])):
            csvrow.append(header[headers][heading].value)
        csvoutput.append(csvrow)


    # current parent_key
    first_col = ''
    for rowIndex in xrange(num_header_rows, sheet.nrows):
        # Read complete row
        row = sheet.row(rowIndex)

        # read name of next parent_key
        next_first_col = row[0].value

        # if empty, use parent_key of previous row
        if next_first_col == u'':
            if first_col == u'':
                raise ValueError('First data row in xls file does not contain a parent_key identification.')
        else:
            first_col = next_first_col

        # concatenate first_col and parent_key
        parent_key = first_col
        for keypart in xrange(1, num_parent_key_cols):
            parent_key = parent_key + parent_key_sep + row[keypart].value

        # check: is this parent_key unique?
        if parent_key in result.keys():
            raise ValueError('Parent key exists multiple times in xls file.')

        # create new dictionary for this parent_key
        result[parent_key] = {}

        # go through columns 2-...
        for colIndex in xrange(num_parent_key_cols, sheet.ncols):
            level = result[parent_key]

            # if we have only 1 header row, there is nothing to do
            final_header = 0
            if num_header_rows > 1:
                # if current header does not yet exist as dictionary entry, then
                # create a new dictionary (note: if there are num_header_rows rows with
                # headings, we have only num_header_rows-1 nested dictionaries)
                for headers in xrange(num_header_rows - 1):
                    # as long as this row is not this row not empty?
                    if not header[headers][colIndex].value == u'':
                        # is the next row empty? if yes, then the current row the final key
                        if header[headers + 1][colIndex].value == u'':
                            final_header = headers
                        # if no, (possibly) create and progress to the next a new nested level
                        else:
                            if not header[headers][colIndex].value in level:
                                level[header[headers][colIndex].value] = {}
                            level = level[header[headers][colIndex].value]

                            # is this the last row? then the next (non-empty) row is the final key
                            final_header = headers + 1


            # save value to lowest nested dictionary
            level[header[final_header][colIndex].value] = row[colIndex].value


        # create row for csv output
        csvrow = [row[i].value for i in  xrange(sheet.ncols)]
        csvoutput.append(csvrow)

    return result, csvoutput