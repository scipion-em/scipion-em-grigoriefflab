# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.em

from .constants import *


_logo = "brandeis_logo.png"


def _getHome(binaryKey, default, paths):
    """ Get the required home path, if not present..
    the default value will be used from EM_ROOT.
    Can join extra paths
    """
    home = os.environ.get('%s_HOME' % binaryKey,
                          os.path.join(os.environ['EM_ROOT'], default))
    return os.path.join(home, *paths)


def _getCtffind4Paths():
    if os.path.exists(_getHome(CTFFIND4, 'ctffind4', ['bin', CTFFIND4_BIN])):
        return ['bin']
    return []


# This is kind of special plugin since the program binaries are distributed
# separated. So, there is not a single "HOME" directory and a single version.
class Plugin(pyworkflow.em.Plugin):
    __programs = {
        CTFFIND: {
            'DEFAULT': 'ctffind',
            '': [CTFFIND_BIN, CTFFINDMP_BIN],  # default program, 1 or 2 exec
            CTFTILT: [CTFTILT_BIN, CTFTILTMP_BIN]
        },
        CTFFIND4: {
            'DEFAULT': 'ctffind4',
            'PATH': _getCtffind4Paths(),  # variable path, depending on 'bin'
            '': [CTFFIND4_BIN],
            'VERSIONS': ['4.0.15', '4.1.5', '4.1.8', V4_1_10]
        },
        FREALIGN: {
            'DEFAULT': 'frealign',
            'PATH': ['bin'],
            '': [FREALIGN_BIN, FREALIGNMP_BIN],
            CALC_OCC: [CALC_OCC_BIN],
            RSAMPLE: [RSAMPLE_BIN],
            'VERSIONS': ['9.07']
        },
        MAGDIST: {
            'DEFAULT': 'mag_distortion',
            'PATH': ['bin'],
            MAGDISTEST: [MAGDISTEST_BIN],
            MAGDISTCORR: [MAGDISTCORR_BIN]
        },
        UNBLUR: {
            'DEFAULT': 'unblur',
            'PATH': ['bin'],
            '': [UNBLUR_BIN],
            'VERSIONS': ['1.0_150529', '1.0.2']
        },
        SUMMOVIE: {
            'DEFAULT': 'summovie',
            'PATH': ['bin'],
            '': [SUMMOVIE_BIN]
        }
    }

    @classmethod
    def __getEntry(cls, binaryKey):
        entry = cls.__programs.get(binaryKey, None)

        if entry is None:
            raise Exception("Binaries for '%S' not found. " % binaryKey)

        return entry

    @classmethod
    def getProgram(cls, binaryKey, programKey='', useMP=False):
        """ Get the whole path of a given program.
        Programs are grouped by binaries.
        If programKey is '', the default one will be taken.
        """
        entry = cls.__getEntry(binaryKey)

        if programKey not in entry:
            raise Exception("Invalid program '%s' for binary %s"
                            % (programKey, binaryKey))

        n = len(entry[programKey])

        if n < 1 or (useMP and n != 2):
            raise Exception("Wrong number of binaries or not support for MP.")

        program = entry[programKey][1 if useMP else 0]
        paths = entry.get('PATH', []) + [program]
        return _getHome(binaryKey, entry['DEFAULT'], paths)

    @classmethod
    def getHome(cls, binaryKey=FREALIGN):
        entry = cls.__getEntry(binaryKey)
        return _getHome(binaryKey, entry['DEFAULT'], entry.get('PATH', []))

    @classmethod
    def getSupportedVersions(cls, binaryKey=FREALIGN):
        """ Return the list of supported binary versions. """
        return cls.__getEntry(binaryKey).get('VERSIONS', [])

    @classmethod
    def getActiveVersion(cls, binaryKey=FREALIGN):
        """ Return the version of the Relion binaries that is currently active.
        In the current implementation it will be inferred from the RELION_HOME
        variable, so it should contain the version number in it. """
        home = cls.getHome(binaryKey)
        for v in cls.getSupportedVersions(binaryKey):
            if v in home:
                return v
        return ''


pyworkflow.em.Domain.registerPlugin(__name__)

# TODO: Remove the following lines when no longer needed
#
# CTFFIND_PATH = join(os.environ[CTFFIND_HOME], CTFFIND3)
# CTFFINDMP_PATH = join(os.environ[CTFFIND_HOME], CTFFIND3MP)
# CTFFIND4_PATH = _getCtffind4()
#
# CTFTILT_PATH = join(os.environ[CTFFIND_HOME], CTFTILT)
# CTFTILTMP_PATH = join(os.environ[CTFFIND_HOME], CTFTILTMP)
#
# FREALIGN_HOME = _getHome(FREALIGN_HOME, 'frealign')
# FREALIGN_PATH = join(FREALIGN_HOME, 'bin', FREALIGN)
# FREALIGNMP_PATH = join(FREALIGN_HOME, 'bin', FREALIGNMP)
#
# CALC_OCC_PATH = join(FREALIGN_HOME, 'bin', CALC_OCC)
# RSAMPLE_PATH = join(FREALIGN_HOME, 'bin', RSAMPLE)
#
# MAGDIST_HOME = _getHome(MAGDIST_HOME, 'mag_distortion')
# MAGDISTEST_PATH = join(MAGDIST_HOME, 'bin', MAGDISTEST_BIN)
# MAGDISTCORR_PATH = join(MAGDIST_HOME, 'bin', MAGDISTCORR_BIN)
#
# UNBLUR_PATH = join(_getHome(UNBLUR_HOME, 'unblur'), 'bin', UNBLUR_BIN)
# SUMMOVIE_PATH = join(_getHome(SUMMOVIE_HOME, 'summovie'), 'bin', SUMMOVIE_BIN)
#
#
# def validateMagDistorsionInstallation():
#     """ Check if the installation of this protocol is correct.
#     Can't rely on package function since this is a "multi package" package
#     Returning an empty list means that the installation is correct
#     and there are not errors. If some errors are found, a list with
#     the error messages will be returned.
#     """
#     missingPaths = []
#
#     if not os.path.exists(MAGDIST_HOME):
#         missingPaths.append("%s : %s" % (MAGDIST_HOME, MAGDIST_HOME))
#     return missingPaths