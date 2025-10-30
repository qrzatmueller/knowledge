#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.10.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/pqrz/Documents/MasterThesis/mt_eq/VerbindungUnt/Simulation/Geometry')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
StahlBlech = geompy.MakeFaceHW(0.15, 0.03, 1)
geompy.TranslateDXDYDZ(StahlBlech, 0.075, 0, 0)
Verbindung = geompy.MakeBoxDXDYDZ(0.02, 0.03, 0.005)
geompy.TranslateDXDYDZ(Verbindung, 0.13, -0.015, 0)
Verbindung_face_33 = geompy.GetSubShape(Verbindung, [33])
Face_1 = geompy.MakeFaceObjHW(Verbindung_face_33, 0.3, 0.03)
geompy.TranslateDXDYDZ(Face_1, 0.14, 0, 0)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( StahlBlech, 'StahlBlech' )
geompy.addToStudy( Verbindung, 'Verbindung' )
geompy.addToStudyInFather( Verbindung, Verbindung_face_33, 'Verbindung:face_33' )
geompy.addToStudy( Face_1, 'Face_1' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
