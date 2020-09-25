#!/usr/bin/env python

## \file p3d2su2_3D.py
#  \brief Python script for converting from plot3D to SU2
#  \author F. Palacios
#  \version 3.2.8 "eagle"
#
# SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

from optparse import OptionParser
import string

# parser=OptionParser()
# parser.add_option("-f", "--file", dest="filename", default="default.p3d",
#                   help="write mesh to FILE", metavar="FILE")
# (options, args)=parser.parse_args()

filename = 'meshes/CST-Test-100.xyz'
# Read the input file
p3d_File = open(filename,"r")

# Read the header
next(p3d_File)
header = p3d_File.readline().replace("\n"," ").replace("\t"," ").split()
print(header)
nNode = int(header[0].strip())
mNode = int(header[1].strip())
lNode = int(header[2].strip())

# Read the body
body = p3d_File.read().replace("\n"," ").replace("\t"," ").split()

p3d_File.close()

# Write the .su2 file
filename = filename.rsplit( ".", 1 )[ 0 ] + ".su2"
su2_File = open(filename,"w")

# Write the header
su2_File.write( "NDIME=3\n" )
su2_File.write( "NELEM=%s\n" % ((lNode-1)*(nNode-1)*(mNode-1)))

# Write the connectivity
iElem = 0
for kNode in range(lNode-1):
    for jNode in range(mNode-1):
        for iNode in range(nNode-1):
            Point0 = kNode*mNode*nNode + jNode*nNode + iNode
            Point1 = kNode*mNode*nNode + jNode*nNode + iNode + 1
            Point2 = kNode*mNode*nNode + (jNode+1)*nNode + (iNode+1)
            Point3 = kNode*mNode*nNode + (jNode+1)*nNode + iNode
            Point4 = (kNode+1)*mNode*nNode + jNode*nNode + iNode
            Point5 = (kNode+1)*mNode*nNode + jNode*nNode + iNode + 1
            Point6 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + (iNode + 1)
            Point7 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + iNode
            su2_File.write( "12 \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n" % (Point0, Point1, Point2, Point3, Point4, Point5, Point6, Point7, iElem) )
            iElem = iElem + 1

# Write the coordinates
nPoint = (nNode)*(mNode)*(lNode)
su2_File.write( "NPOIN=%s\n" % ((nNode)*(mNode)*(lNode)))
iPoint = 0
for kNode in range(lNode):
    for jNode in range(mNode):
        for iNode in range(nNode):
            XCoord = body[kNode*(mNode*nNode) + jNode*nNode + iNode]
            YCoord = body[(nNode*mNode*lNode) + kNode*(mNode*nNode) + jNode*nNode + iNode]
            ZCoord = body[2*(nNode*mNode*lNode) + kNode*(mNode*nNode) + jNode*nNode + iNode]
            su2_File.write( "%s \t %s \t %s \t %s\n" % (XCoord, YCoord, ZCoord, iPoint) )
            iPoint = iPoint + 1

# Write the boundaries
su2_File.write( "NMARK=6\n" )

su2_File.write( "MARKER_TAG= Xplane_0\n" )
elem = (mNode-1)*(lNode-1)
su2_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-1):
    for kNode in range(lNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % (jNode*nNode + (nNode - 1) + kNode*nNode*mNode,  (jNode + 1)*nNode + (nNode - 1) + kNode*nNode*mNode,  (jNode + 1)*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode, jNode*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode ) )

su2_File.write( "MARKER_TAG= Xplane_1\n" )
elem = (mNode-1)*(lNode-1)
su2_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-2, -1, -1):
    for kNode in range(lNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % ((jNode + 1)*nNode + kNode*nNode*mNode, jNode*nNode + kNode*nNode*mNode, jNode*nNode+ (kNode+1)*nNode*mNode, (jNode + 1)*nNode+ (kNode+1)*nNode*mNode ) )

su2_File.write( "MARKER_TAG= Yplane_0\n" )
elem = (nNode-1)*(lNode-1)
su2_File.write( "MARKER_ELEMS=%s\n" % elem )
for iNode in range(nNode-1):
    for kNode in range(lNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % (iNode + kNode*nNode*mNode, iNode + (kNode+1)*nNode*mNode, iNode + 1 + (kNode+1)*nNode*mNode, iNode + 1 + kNode*nNode*mNode) )

su2_File.write( "MARKER_TAG= Yplane_1\n" )
elem = (nNode-1)*(lNode-1)
su2_File.write( "MARKER_ELEMS=%s\n" % elem )
for iNode in range(nNode-1):
    for kNode in range(lNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % ((nNode*mNode - 1) - iNode + kNode*nNode*mNode,  (nNode*mNode - 1) - iNode + (kNode+1)*nNode*mNode, (nNode*mNode - 1) - (iNode + 1) + (kNode+1)*nNode*mNode, (nNode*mNode - 1) - (iNode + 1) + kNode*nNode*mNode) )

su2_File.write( "MARKER_TAG= Zplane_0\n" )
elem = (nNode-1)*(mNode-1);
su2_File.write( "MARKER_ELEMS=%s\n" % elem)
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % (jNode*nNode + iNode, jNode*nNode + (iNode+1), (jNode + 1)*nNode + (iNode + 1), (jNode + 1)*nNode + iNode) )

su2_File.write( "MARKER_TAG= Zplane_1\n" )
elem = (nNode-1)*(mNode-1)
su2_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        su2_File.write( "9 \t %s \t %s \t %s \t %s\n" % (nNode*mNode*(lNode - 1) + jNode*nNode + iNode, nNode*mNode*(lNode - 1) + jNode*nNode + iNode + 1, nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + (iNode + 1), nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + iNode) )

su2_File.close()