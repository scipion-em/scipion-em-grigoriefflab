# **************************************************************************
# *
# * Authors:     Josue Gomez BLanco (jgomez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol for CTF estimation with ctffind3
"""


from pyworkflow.utils.path import replaceBaseExt, cleanPattern
from pyworkflow.em import *
from pyworkflow.em.packages.grigoriefflab import CTFFIND_PATH, CTFFIND4_PATH
from pyworkflow.em.packages.grigoriefflab.convert import readCtfModel

class ProtCTFFind(ProtCTFMicrographs):
    """Estimates CTF on a set of micrographs
    using either ctffind3 or ctffind4 program"""
    _label = 'ctffind'
    
    
    def _defineProcessParams(self, form):
        form.addParam('useCftfind4', BooleanParam, default=False,
              label="use ctffind4 to estimate the CTF?",
              help='If is true, the protocol will use ctffind4 instead of ctffind3')
        form.addParam('astigmatism', FloatParam, default=100.0,
              label='Expected (tolerated) astigmatism', expertLevel=LEVEL_ADVANCED,
              condition='useCftfind4', )
        form.addParam('findPhaseShift', BooleanParam, default=False,
              label="Find additional phase shift?", condition='useCftfind4',
              expertLevel=LEVEL_ADVANCED,)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _estimateCTF(self, micFn, micDir):
        """ Run ctffind, 3 or 4, with required parameters """
        # Create micrograph dir 
        makePath(micDir)
        downFactor = self.ctfDownFactor.get()
        
        micFnMrc = self._getTmpPath(replaceBaseExt(micFn, 'mrc'))
        if downFactor != 1:
            #Replace extension by 'mrc' cause there are some formats that cannot be written (such as dm3)
            self.runJob("xmipp_transform_downsample","-i %s -o %s --step %f --method fourier" % (micFn, micFnMrc, downFactor))
            self._params['scannedPixelSize'] = self.inputMicrographs.get().getScannedPixelSize() * downFactor
        else:
            micFnMrc = self._getTmpPath(replaceBaseExt(micFn, "mrc"))
            ImageHandler().convert(micFn, micFnMrc, DT_FLOAT)
        
        # Update _params dictionary
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = self._getCtfOutPath(micDir)
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception:
            raise
        cleanPath(micFnMrc)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def _restimateCTF(self, id):
        """ Run ctffind3 with required parameters """

        ctfModel = self.recalculateSet[id]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)
        
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        
        cleanPath(out)
        cleanPath(psdFile)
        micFnMrc = self._getTmpPath(replaceBaseExt(micFn, "mrc"))
        ImageHandler().convert(micFn, micFnMrc, DT_FLOAT)

        # Update _params dictionary
        self._prepareRecalCommand(ctfModel)
        self._params['micFn'] = micFnMrc
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = out
        self._params['ctffindPSD'] = psdFile
        
        try:
            self.runJob(self._program, self._args % self._params)
        except Exception:
            raise
        cleanPattern(micFnMrc)
    
    def _createNewCtfModel(self, mic):
        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)
        ctfModel2 = CTFModel()
        
        if not self.useCftfind4:
            readCtfModel(ctfModel2, out)
        else:
            readCtfModel(ctfModel2, out, True)
        
        ctfModel2.setPsdFile(psdFile)
        ctfModel2.setMicrograph(mic)
        return ctfModel2
    
    def _createOutputStep(self):
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(self.inputMics)
        defocusList = []
        
        for fn, micDir, mic in self._iterMicrographs():
            samplingRate = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(samplingRate)
            psdFile = self._getPsdPath(micDir)
            out = self._getCtfOutPath(micDir)
            
            ctfModel = CTFModel()
        
            if not self.useCftfind4:
                readCtfModel(ctfModel, out)
            else:
                readCtfModel(ctfModel, out, True)
            
            ctfModel.setPsdFile(psdFile)
            ctfModel.setMicrograph(mic)
            
            defocusList.append(ctfModel.getDefocusU())
            defocusList.append(ctfModel.getDefocusV())
            ctfSet.append(ctfModel)
        
        self._defocusMaxMin(defocusList)
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(self.inputMics, ctfSet)
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        ctffind = join(os.environ['CTFFIND_HOME'], 'ctffind3.exe')
        if not exists(ctffind):
            errors.append('Missing ctffind3.exe')
        return errors
    
    def _citations(self):
        return ['Mindell2003']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using CtfFind [Midell2003]. " % self.getObjectTag('inputMicrographs')
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')
        
        return [methods]
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _prepareCommand(self):
        sampling = self.inputMics.getSamplingRate() * self.ctfDownFactor.get()
        # Convert digital frequencies to spatial frequencies
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / self._params['lowRes']
        if self._params['lowRes'] > 50:
            self._params['lowRes'] = 50
        self._params['highRes'] = sampling / self._params['highRes']
        self._params['step_focus'] = 500.0
        if not self.useCftfind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
            else:
                self._params['phaseShift'] = "no"
            self._argsCtffind4()
    
    def _prepareRecalCommand(self, ctfModel):
        line = ctfModel.getObjComment().split()
        self._defineRecalValues(ctfModel)
        # get the size and the image of psd
    
        imgPsd = ctfModel.getPsdFile()
        imgh = ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)
        
        mic = ctfModel.getMicrograph()
        micDir = self._getMicrographDir(mic)
        
        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        self._params['step_focus'] = 1000.0
        self._params['sampling'] = sampling
        self._params['lowRes'] = sampling / float(line[3])
        self._params['highRes'] = sampling / float(line[4])
        self._params['minDefocus'] = min([float(line[0]), float(line[1])])
        self._params['maxDefocus'] = max([float(line[0]), float(line[1])])
        self._params['windowSize'] = size
        
        if not self.useCftfind4:
            self._argsCtffind3()
        else:
            self._params['astigmatism'] = self.astigmatism.get()
            if self.findPhaseShift:
                self._params['phaseShift'] = "yes"
            else:
                self._params['phaseShift'] = "no"
            self._argsCtffind4()
    
    def _argsCtffind3(self):
        self._program = 'export NATIVEMTZ=kk ; ' + CTFFIND_PATH
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
eof
"""
    
    def _argsCtffind4(self):
            
        self._program = CTFFIND4_PATH
        self._args = """ << eof
%(micFn)s
%(ctffindPSD)s
%(sampling)f
%(voltage)f
%(sphericalAberration)f
%(ampContrast)f
%(windowSize)d
%(lowRes)f
%(highRes)f
%(minDefocus)f
%(maxDefocus)f
%(step_focus)f
%(astigmatism)f
%(phaseShift)s
eof
"""            
    
    def _getPsdPath(self, micDir):
        return join(micDir, 'ctfEstimation.mrc')
    
    def _getCtfOutPath(self, micDir):
        return join(micDir, 'ctfEstimation.txt')
    
#     def _parseOutput(self, filename):
#         """ Try to find the output estimation parameters
#         from filename. It search for a line containing: Final Values.
#         """
#         if not self.useCftfind4:
#             return parseCtffindOutput(filename)
#         else:
#             return parseCtffind4Output(filename)
    
#     def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
#         ctf = CTFModel()
#         ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
#         ctf.setPsdFile(psdFile)
#         
#         return ctf
