# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import os
import sys
import pyworkflow.utils as pwutils
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_2

from grigoriefflab import Plugin
from grigoriefflab.constants import (CTFFIND, CTFTILT)
from grigoriefflab.convert import readCtfModel, parseCtftiltOutput


class ProtCTFTilt(em.ProtCTFMicrographs):
    """
    Estimates CTF on a set of tilted micrographs
    using ctftilt program.
    
    """
    _label = 'ctftilt'
    _lastUpdateVersion = VERSION_1_2

    @classmethod
    def validateInstallation(cls):
        """ Check if the installation of this protocol is correct.
        Can't rely on package function since this is a "multi package" package
        Returning an empty list means that the installation is correct
        and there are no errors. If some errors are found, a list with
        the error messages will be returned.
        """
        missingPaths = []

        # FIXME

        if not os.path.exists(Plugin.getHome(CTFTILT)):
            missingPaths.append("%s : ctffind3/ctftilt installation not found"
                                     % CTFTILT)

        # ctftilt = cls._getProgram()
        #
        # if not os.path.exists(ctftilt):
        #     missingPaths.append("%s : ctffind3/ctftilt installation not found"
        #                         % ctftilt)
        return missingPaths

    def _defineProcessParams(self, form):
        form.addParam('astigmatism', params.FloatParam, default=100.0,
                      label='Expected astigmatism (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Expected amount of astigmatism in Angstrom. ')

        line = form.addLine('Tilt angle',
                            help='Expected tilt angle value and its '
                                 'uncertainty in degrees.')
        line.addParam('tiltA', params.FloatParam, default=0.,
                      label='Expected value')
        line.addParam('tiltR', params.FloatParam, default=5.,
                      label='Uncertainty')

    # --------------------------- STEPS functions -----------------------------
    def _estimateCTF(self, mic, *args):
        """ Run ctftilt with required parameters """
        micFn = mic.getFileName()
        micDir = self._getTmpPath('mic_%04d' % mic.getObjId())
        doneFile = os.path.join(micDir, 'done.txt')

        if self.isContinued() and os.path.exists(doneFile):
            return

        try:
            # Create micrograph dir
            pwutils.makePath(micDir)
            downFactor = self.ctfDownFactor.get()
            scannedPixelSize = self.inputMicrographs.get().getScannedPixelSize()
            micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, 'mrc'))

            if downFactor != 1:
                # Replace extension by 'mrc' because there are some formats
                # that cannot be written (such as dm3)
                em.ImageHandler().scaleFourier(micFn, micFnMrc, downFactor)
                self._params['scannedPixelSize'] = scannedPixelSize * downFactor
            else:
                ih = em.ImageHandler()
                if ih.existsLocation(micFn):
                    micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
                    ih.convert(micFn, micFnMrc, em.DT_FLOAT)
                else:
                    print >> sys.stderr, "Missing input micrograph %s" % micFn

        except Exception as ex:
            print >> sys.stderr, "Some error happened: %s" % ex
            import traceback
            traceback.print_exc()

        try:
            program, args = self._getCommand(micFn=micFnMrc,
                                             ctftiltOut=self._getCtfOutPath(micDir),
                                             ctftiltPSD=self._getPsdPath(micDir))
            self.runJob(program, args)
        except Exception as ex:
            print >> sys.stderr, "ctftilt has failed with micrograph %s" % micFnMrc

        # Let's notify that this micrograph have been processed
        # just creating an empty file at the end (after success or failure)
        open(doneFile, 'w')
        # Let's clean the temporary mrc micrographs
        pwutils.cleanPath(micFnMrc)

    def _restimateCTF(self, ctfId):
        """ Run ctftilt with required parameters """

        ctfModel = self.recalculateSet[ctfId]
        mic = ctfModel.getMicrograph()
        micFn = mic.getFileName()
        micDir = self._getMicrographDir(mic)

        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        pwutils.cleanPath(out)
        micFnMrc = self._getTmpPath(pwutils.replaceBaseExt(micFn, "mrc"))
        em.ImageHandler().convert(micFn, micFnMrc, em.DT_FLOAT)
        pwutils.cleanPath(psdFile)
        try:
            program, args = self._getRecalCommand(
                ctfModel, micFn=micFnMrc, ctftiltOut=out, ctftiltPSD=psdFile)
            self.runJob(program, args)
        except Exception as ex:
            print >> sys.stderr, "ctftilt has failed with micrograph %s" % micFnMrc
        pwutils.cleanPattern(micFnMrc)

    def _createCtfModel(self, mic, updateSampling=True):
        #  When downsample option is used, we need to update the
        # sampling rate of the micrograph associated with the CTF
        # since it could be downsampled
        if updateSampling:
            newSampling = mic.getSamplingRate() * self.ctfDownFactor.get()
            mic.setSamplingRate(newSampling)

        micDir = self._getMicrographDir(mic)
        out = self._getCtfOutPath(micDir)
        psdFile = self._getPsdPath(micDir)

        ctfModel = em.CTFModel()
        readCtfModel(ctfModel, out, ctf4=False, ctfTilt=True)
        ctfModel.setPsdFile(psdFile)
        ctfModel.setMicrograph(mic)

        return ctfModel

    def _createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        ctftilt = self._getProgram()

        if not os.path.exists(ctftilt):
            errors.append('Missing %s' % ctftilt)

        return errors

    def _citations(self):
        return ['Mindell2003']

    def _methods(self):
        if self.inputMicrographs.get() is None:
            return ['Input micrographs not available yet.']
        methods = "We calculated the CTF of %s using CTFTilt. " % self.getObjectTag('inputMicrographs')
        methods += self.methodsVar.get('')
        methods += 'Output CTFs: %s' % self.getObjectTag('outputCTF')

        return [methods]

    # -------------------------- UTILS functions ------------------------------
    def _getProgram(self):
         return Plugin.getProgram(CTFTILT,
                                  useMP=self.numberOfThreads > 1)

    def _getCommand(self, **kwargs):
        params = self.getCtfParamsDict()
        # Convert digital frequencies to spatial frequencies
        sampling = params['samplingRate']
        params['lowRes'] = sampling / params['lowRes']
        if params['lowRes'] > 50:
            params['lowRes'] = 50
        params['highRes'] = sampling / params['highRes']
        params['astigmatism'] = self.astigmatism.get()
        params['step_focus'] = 500.0
        params['pixelAvg'] = 1  # set to 1 since we have our own downsampling
        params['tiltAngle'] = self.tiltA.get()
        params['tiltR'] = self.tiltR.get()
        params.update(kwargs)
        return self._getCommandFromParams(params)

    def _getRecalCommand(self, ctfModel, **kwargs):
        line = ctfModel.getObjComment().split()
        params = self.getRecalCtfParamsDict(ctfModel)
        # get the size and the image of psd

        imgPsd = ctfModel.getPsdFile()
        imgh = em.ImageHandler()
        size, _, _, _ = imgh.getDimensions(imgPsd)

        mic = ctfModel.getMicrograph()

        # Convert digital frequencies to spatial frequencies
        sampling = mic.getSamplingRate()
        params['step_focus'] = 1000.0
        params['sampling'] = sampling
        params['lowRes'] = sampling / float(line[3])
        params['highRes'] = sampling / float(line[4])
        params['minDefocus'] = min([float(line[0]), float(line[1])])
        params['maxDefocus'] = max([float(line[0]), float(line[1])])
        params['astigmatism'] = self.astigmatism.get()
        params['windowSize'] = size
        params['pixelAvg'] = 1  # set to 1 since we have our own downsampling
        params['tiltAngle'] = self.tiltA.get()
        params['tiltR'] = self.tiltR.get()
        params.update(kwargs)

        return self._getCommandFromParams(params)

    def _useThreads(self):
        return self.numberOfThreads > 1

    def _getCommandFromParams(self, params):
        program = 'export NATIVEMTZ=kk ; '
        if self._useThreads():
            program += 'export NCPUS=%d ;' % self.numberOfThreads
        program += self._getProgram()
        args = """   << eof > %(ctftiltOut)s
%(micFn)s
%(ctftiltPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f,%(pixelAvg)d
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f,%(astigmatism)f,%(tiltAngle)f,%(tiltR)f
eof
"""
        return program, args % params

    def _getPsdPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.mrc')

    def _getCtfOutPath(self, micDir):
        return os.path.join(micDir, 'ctfEstimation.txt')

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It searches for a line containing: Final Values.
        """
        return parseCtftiltOutput(filename)

    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = em.CTFModel()
        ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
        ctf.setPsdFile(psdFile)

        return ctf

    def _summary(self):
        summary = em.ProtCTFMicrographs._summary(self)
        if hasattr(self, 'outputCTF'):
            ctfs = self.outputCTF
            for ctf in ctfs:
                angle = float(ctf._ctftilt_tiltAngle)
                axis = float(ctf._ctftilt_tiltAxis)
                summary.append('Estimated tilt parameters:\n - tilt angle _%0.2f_\n'
                               ' - tilt axis _%0.2f_' %
                               (angle, axis))
                summary.append('If you think that tilt angle should have an '
                               'opposite sign than reported, use the following '
                               'values:\n - tilt angle _%0.2f_\n'
                               ' - tilt axis _%0.2f_' %
                               (-angle, axis + 180.0) )

        return summary
