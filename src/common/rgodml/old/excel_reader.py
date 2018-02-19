'''
Created on Oct 30, 2013

@author: denker
'''
import re
import xlrd
import datetime

def xls_sheet_reader(filename, sheetname, num_parent_key_cols, parent_key_sep, num_header_rows):
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
    xlsbook = xlrd.open_workbook(filename)
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

def datefloat2date(datefloat):
    """
    Converts date float numbers from Excel to python datetimes.

    Args:
        datefloat (float):
            Float number representing a certain date.

    Returns:
        date (datetime.date)
    """
    y1900 = datetime.datetime(1900,1,1)
    y1970 = datetime.datetime(1970,1,1)

    from1900to1970 = y1970 - y1900 + datetime.timedelta(days=2)

    date = datetime.date.fromtimestamp(int(datefloat) * 86400) - from1900to1970

    return date

def string2list(string, dtype):
    """
    Converts string from Excel to a python list with values of type int or float.

    Args:
        string (string):
            String representing a list of integers or floats.
        dtype:
            Data type of values in wanted list.

    Returns:
        list (of int or float)
    """
    if dtype == 'int':
        vlist = [int(i) for i in re.findall(r'\d+', string)]
    elif dtype == 'float':
        vlist = [float(i) for i in re.findall(r'\d.\d+', string)]
    else:
        raise IOError("The data type can not be converted into a list. "+
                      "Possible data types are 'int' or 'float'.")

    return vlist

#if __name__ == '__main__':

#    filename = '/home/zehl/datasets/marseille/DataGrasp/MetaDataLilou/source/SubsessionInfo/Lilou_spikesortinginfo'
#    dictout, csvout = xls_sheet_reader(filename=filename + '.xls', sheetname="Sheet1", num_parent_key_cols=2, parent_key_sep='-', num_header_rows=2)
#    print dictout
#    import csv
#    csvbook = csv.writer(open(filename + '.csv', 'wb'), delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
#    csvbook.writerows(csvout)

#    filename = '/home/zehl/datasets/marseille/DataGrasp/MetaDataLilou/source/SubsessionInfo/Lilou_subsessioninfo'
#    dictout, csvout = xls_sheet_reader(filename=filename + '.xls', sheetname="Sheet1", num_parent_key_cols=2, parent_key_sep='-', num_header_rows=5)
#    print dictout
#    import csv
#    csvbook = csv.writer(open(filename + '.csv', 'wb'), delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
#    csvbook.writerows(csvout)

#    filename = '/home/zehl/datasets/marseille/DataGrasp/MetaDataGeneral/GeneralInfo'
#    dictout, csvout = xls_sheet_reader(filename=filename + '.xls', sheetname="Sheet1", num_parent_key_cols=2, parent_key_sep='', num_header_rows=1)
#    print dictout
#    import csv
#    csvbook = csv.writer(open(filename + '.csv', 'wb'), delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
#    csvbook.writerows(csvout)


