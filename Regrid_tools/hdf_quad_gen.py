#!/usr/bin/env python

# Python code that will generate an appropriate hdf_quadrangle file given a text file with the field names to grid

import os.path
import sys

fields_file = '/Users/Josh/Documents/MATLAB/BEHR/OMI_SP/OMI_fields.txt'
cellfields_file = '/Users/Josh/Documents/MATLAB/BEHR/OMI_SP/OMI_cellfields.txt'
hdf_quad_template = 'hdf_quadrangle_template.m'
hdf_quad_save = '/Users/Josh/Documents/MATLAB/BEHR/OMI_SP/hdf_quadrangle_OMI.m'

# Check if the save file exists, if the user declines to overwrite, abort
if os.path.isfile(hdf_quad_save):
    ans = raw_input('Save file exists. Overwrite? y/n [n]  ')
    if (ans != 'y') & (ans != 'Y'):
        print('Aborting script.  Change value of hdf_quad_save to avoid overwriting existing file.')
        sys.exit()

# Save all field names to the variable fields
with open(fields_file,'r') as f:
    fields = []
    for line in f:
        fields.append(line.strip())
        
with open(cellfields_file,'r') as f2:
	cellfields = []
	for line in f2:
		cellfields.append(line.strip())
        
numfields = len(fields)
numcellfields = len(cellfields)

# Get the function name out of the file name
func_name = os.path.basename(hdf_quad_save).split('.')[0]

# Open the template file and the file to be saved
f_template = open(hdf_quad_template,'r')
f_save = open(hdf_quad_save,'w')

# Copy each line from the template file, unless the first two non-whitespace characters are %$
# The % sign signifies a Matlab comment, the $ is not generally used in Matlab, but here indicates
#   that this line represents a pattern which the field names are to be inserted into.
#   Field names will then replace any instance of $field in the rest of the line.  For example:
#
# As an example, imagine we wanted to grid 3 fields: NO2, Latitude, and Longitude.  In the Matlab script:
#
#       %$ $field_val = Data.$field;
#
#   would become:
#
#       NO2_val = Data.NO2;
#       Latitude_val = Data.Latitude;
#       Longitude_val = Data.Longitude;
#
# The wildcard $keyfield will always be replaced with the first field in the list.  Given the fields above, the line:
#
#       %$ if ~isnan(Data.$keyfield) && ~isnan(Data.$field)
#
#   would become
#
#       if ~isnan(Data.NO2) && ~isnan(Data.NO2)
#       if ~isnan(Data.NO2) && ~isnan(Data.Latitude)
#       if ~isnan(Data.NO2) && ~isnan(Data.Longitude)
#
# A line whose first non-whitespace characters are %# will only have $keyfield replaced and will not be replicated, i.e.
#
#       %# if ~isnan(Data.$keyfield)
#
#   becomes
#
#       if ~isnan(Data.NO2)
#

for line in f_template:
    test_line = line.lstrip()
    line = line.replace('hdf_quadrangle_template',func_name)
    if test_line[0:3] == '%$f':
        newline = line.replace('%$f ','')
        for a in range(numfields):
            lineout = newline.replace('$field',fields[a])
            lineout = lineout.replace('$keyfield',fields[0])
            f_save.write(lineout)
    elif test_line[0:3] == '%$c':
    	newline = line.replace('%$c ','')
    	for a in range(numcellfields):
    		lineout = newline.replace('$cellfield',cellfields[a])
    		lineout = lineout.replace('$keyfield',fields[0])
    		f_save.write(lineout)
    elif test_line[0:2] == '%#':
        lineout = line.replace('%#','')
        lineout = lineout.replace('$keyfield',fields[0])
        f_save.write(lineout)
    else:
        f_save.write(line)
        
f_template.close()
f_save.close()
