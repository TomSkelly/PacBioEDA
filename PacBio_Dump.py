#!/usr/bin/env python

# Copyright (C) 2011 Genome Research Limited -- See full notice at end
# of module.

# Python script to dump a generic hdf5 file. Doesn't know anything
# about PacBio formats. Run it using run_py. First (only) parameter is
# input hdf5 file. Output is to stdout, redirect with >.

import sys
import h5py

def print_attrs (obj, depth):

    for attr in obj.attrs.items():

        print '    ' * (depth+1), '* ', attr[0], '=', attr[1]

def query (group, depth):

    saved_groups = []

    for item in group.items():

        if isinstance (item[1], h5py.highlevel.Dataset):
        
            print '    ' * depth, 'Dataset: ', item[0], item[1].dtype, item[1].shape
            print_attrs (item[1], depth)
            
        elif isinstance (item[1], h5py.highlevel.Group):

            saved_groups.append (item)           # print groups after datasets

        else:

            print '    ' * depth, 'Unknown: ', item[0], item[1].__class__

    for saved_item in saved_groups:

        print '    ' * depth, 'Group:   ', saved_item[0]
        print_attrs (saved_item[1], depth)
            
        query (saved_item[1], depth+1)


def main ():

    infile_name = sys.argv[1]

    print "\nFile: ", infile_name, "\n"

    infile = h5py.File (infile_name, 'r')


####    top = h5py.Group (infile, '/')
    top = infile

    query (top, 0)

    infile.close()

main()


# Copyright (C) 2011 Genome Research Limited
#
# This library is free software. You can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
