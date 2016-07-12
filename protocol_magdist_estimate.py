# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
# *              Grigory Sharov (sharov@igbmc.fr)
# *
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

#from os.path import exists, join

import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils.path as pwutils
from pyworkflow.em.protocol import ProtPreprocessMicrographs

from grigoriefflab import MAGDISTEST_PATH
from convert import parseMagEstOutput


class ProtMagDistEst(ProtPreprocessMicrographs):
    """ mag_distortion_estimate automatically estimates anisotropic magnification
    distortion from a set of images of a standard gold shadowed diffraction grating
    """    
    _label = 'mag_distortion estimate'

    def __init__(self, **args):
        ProtPreprocessMicrographs.__init__(self, **args)
        self.stepsExecutionMode = cons.STEPS_PARALLEL

    # --------------------------- DEFINE params functions --------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Preprocess')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label="Input micrographs", important=True,
                      help='Select the SetOfMicrograph containing ~20 images of different areas of polycrystalline gold')

        form.addParam('scaleFactor', params.FloatParam, default=0.03, expertLevel=params.LEVEL_ADVANCED,
                      label='Scale factor',
                      help='Maximum allowed scale factor.')
        form.addParam('scaleStep', params.FloatParam, default=0.0005, expertLevel=params.LEVEL_ADVANCED,
                      label='Scale step',
                      help='Step size for the scale search.')

        line = form.addLine('Resolution limit', expertLevel=params.LEVEL_ADVANCED,
                            help='Resolution limits for the search.')
        line.addParam('lowRes', params.FloatParam, default=2.5, label='Low')
        line.addParam('highRes', params.FloatParam, default=2.1, label='High')

        line = form.addLine('Angle range (deg)', expertLevel=params.LEVEL_ADVANCED,
                            help='Allowed angle range for the search.')
        line.addParam('minAng', params.FloatParam, default=0.0, label='Min')
        line.addParam('maxAng', params.FloatParam, default=180.0, label='Max')

        form.addParam('angStep', params.FloatParam, default=0.1, expertLevel=params.LEVEL_ADVANCED,
                      label='Angular step (deg)',
                      help='Step size for the angle search.')

        line = form.addLine('Filter radius (freq.)', expertLevel=params.LEVEL_ADVANCED,
                            help='Filter radius for the amplitude bandpass filter.')

        line.addParam('lowp', params.FloatParam, default=0.2, label='Low-pass')
        line.addParam('highp', params.FloatParam, default=0.01, label='High-pass')

        form.addParam('box', params.IntParam, default=512, expertLevel=params.LEVEL_ADVANCED,
                      label='Amplitude box size',
                      help='Box size for the calculated amplitudes.')


        form.addParallelSection(threads=2, mpi=0)

    def _defineInputs(self):
        """ Store some of the input parameters in a dictionary for
        an easy replacement in the programs command line.
        """
        self.params = {'scaleFactor': self.scaleFactor.get(),
                       'scaleStep': self.scaleStep.get(),
                       'lowRes': self.lowRes.get(),
                       'highRes': self.highRes.get(),
                       'minAng': self.minAng.get(),
                       'maxAng': self.maxAng.get(),
                       'angStep': self.angStep.get(),
                       'lowp': self.lowp.get(),
                       'highp': self.highp.get(),
                       'box': self.box.get(),
                       'pixSize': self.inputMicrographs.get().getSamplingRate(),
                       'nthr': self.numberOfThreads.get()}

    # --------------------------- INSERT steps functions --------------------------------------------

    def _prepareStack(self):
        """ Convert input micrographs into a single mrc stack """
        inputMics = self.inputMicrographs.get()
        stackFn = self._getExtraPath('input_stack.mrcs')
        stackFnMrc = self._getExtraPath('input_stack.mrc')

        for fn in inputMics:
            cmd = ' %s %s' %(fn, stackFn)
            self.runJob('e2proc2d.py', cmd)

        pwutils.moveFile(stackFn, stackFnMrc)


    def _processStack(self):
        self._prepareStack()
        stackFnMrc = self._getExtraPath('input_stack.mrc')
        spectraFn = self._getExtraPath('output_amp.mrc')
        rotAvgFn = self._getExtraPath('output_amp_rot.mrc')
        spectraCorrFn = self._getExtraPath('output_amp_corrected.mrc')
        logFn = self._getExtraPath('mag_dist_estimation.log')

        self.params = {'stackFnMrc': stackFnMrc,
                       'spectraFn': spectraFn,
                       'rotAvgFn': rotAvgFn,
                       'spectraCorrFn': spectraCorrFn,
                       'logFn': logFn}

        self._argsMagDist()

        try:
            self.runJob(self._program, self._args % params)

        except:
            print("ERROR: Mag. distortion estimation for %s failed\n" % stackFnMrc)

    # --------------------------- STEPS functions ---------------------------------------------------

    def createOutputStep(self):
        pass
        #self._defineOutputs(outputMicrographs=None)
        #self._defineTransformRelation(self.inputMics, outputMics)

    # --------------------------- INFO functions ----------------------------------------------------

    def _validate(self):
        validateMsgs = []

        return validateMsgs

    def _citations(self):
        return ["Grant2015"]

    def _summary(self):
        summary = []
        for value in self._parseOutputLog():
            summary.append('dist_amount, dist_angle, major_axis, minor_axis= ' % value)

        return summary

    def _methods(self):
        txt = "The micrographs in set %s have " % self.getObjectTag('inputMicrographs')

        return txt

    # --------------------------- UTILS functions --------------------------------------------
    def _parseOutputLog(self):
        """ Return the distortion amount, angle and two scale params. """
        fnOut = self._getExtraPath('mag_dist_estimation.log')

        return parseMagEstOutput(fnOut)

    def _argsMagDist(self):
        self._program = 'export NCPUS=%(nthr)d ; ' + MAGDISTEST_PATH
        self._args = """   << eof > %(logFn)s
%(stackFnMrc)s
%(spectraFn)s
%(rotAvgFn)s
%(spectraCorrFn)s
%(pixSize)f
YES
%(lowRes)f
%(highRes)f
%(scaleFactor)f
%(scaleStep)f
%(minAng)f
%(maxAng)f
%(angStep)f
%(lowp)f
%(highp)f
%(box)d
eof
"""
