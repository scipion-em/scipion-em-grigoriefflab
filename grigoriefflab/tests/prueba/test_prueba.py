# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
try:
 from itertools import izip
except:
    izip = zip

from pwem.objects import SetOfParticles
from pyworkflow.tests import *
from pwem.protocols import (ProtImportParticles, ProtImportVolumes, ProtImportMicrographs)

from grigoriefflab import *
from grigoriefflab.protocols import *


class TestBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')
        cls.volFn = cls.dataset.getFile('vol2')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        kwargs = {
            'filesPath': pattern,
            'magnification': magnification,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration
        }

        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport = ProtImportMicrographs(**kwargs)
        cls.launchProtocol(cls.protImport)

        # Check that input micrographs have been imported
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)

        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate,
                         importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                        )
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        if importFrom == ProtImportParticles.IMPORT_FROM_SCIPION:
            objLabel = 'from scipion (particles)'
        elif importFrom == ProtImportParticles.IMPORT_FROM_FILES:
            objLabel = 'from file (particles)'

        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel=objLabel,
                                         filesPath=pattern,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         importFrom=importFrom)

        cls.launchProtocol(cls.protImport)
        # Check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. '
                            'outputParticles is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=1.237,
                                       voltage=300,
                                       sphericalAberration=2,
                                       scannedPixelSize=None,
                                       magnification=56000)

    @classmethod
    def runImportMicrographRCT(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=2.28,
                                       voltage=100,
                                       sphericalAberration=2.9,
                                       scannedPixelSize=None,
                                       magnification=50000)

    @classmethod
    def runImportParticleGrigorieff(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportParticles(pattern,
                                      samplingRate=4.,
                                      checkStack=True,
                            importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)
    @classmethod
    def runImportVolumesGrigorieff(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportVolumes(pattern,
                                    samplingRate=4.,
                                    importFrom=ProtImportParticles.IMPORT_FROM_FILES)


class TestImportParticlesYunior(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('grigorieff')

    def test_import(self):
        parFile = self.dataset.getFile('particles/particles_iter_002.par')
        stackFile = self.dataset.getFile('particles/particles.mrc')

        protImport = self.newProtocol(ProtImportParticles,
                                      objLabel='import parfile & stk',
                                      parFile=parFile,
                                      stackFile=stackFile,
                                      samplingRate=9.90,
                                      importFrom=ProtImportParticles.IMPORT_FROM_FREALIGN)

        self.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        outputParticles = getattr(protImport, 'outputParticles', None)

        if outputParticles is None:
            raise Exception('Import failed. Par file: %s' % parFile)

        self.assertTrue(outputParticles.getSize() == 180)

        goldFile = self.dataset.getFile('particles/particles.sqlite')
        goldSet = SetOfParticles(filename=goldFile)

        for p1, p2 in izip(goldSet,
                           outputParticles.iterItems(orderBy=['_micId', 'id'],
                                                     direction='ASC')):
            m1 = p1.getTransform().getMatrix()
            m2 = p2.getTransform().getMatrix()
            self.assertTrue(np.allclose(m1, m2, atol=0.01))

        self.assertTrue(outputParticles.hasCTF())
        self.assertTrue(outputParticles.hasAlignmentProj())