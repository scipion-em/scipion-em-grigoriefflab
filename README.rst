
============================
Grigoriefflab scipion plugin
============================
Plugin to use Grigorieff Lab (not CisTEM) programs within the Scipion framework.

=====
Setup
=====

- **Install this plugin:**

.. code-block::

    scipion installp -p scipion-em-grigoriefflab

Alternatively, in devel mode:

.. code-block::

    scipion installp -p local/path/to/scipion-em-grigoriefflab --devel

.. image:: http://scipion-test.cnb.csic.es:9980/badges/grigoriefflab_devel.svg

=====
Tests
=====
.. code-block::

 scipion test grigoriefflab.tests.test_protocols_grigoriefflab_movies
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab_movies.TestUnblur
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab_movies.TestSummovie
 scipion test grigoriefflab.tests.test_protocols_grigoriefflab_magdist
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab_magdist.TestMagDist
 scipion test grigoriefflab.tests.test_protocols_grigoriefflab
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab.TestImportParticles
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab.TestFrealignRefine
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab.TestFrealignClassify
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab.TestCtftilt
   scipion test grigoriefflab.tests.test_protocols_grigoriefflab.TestCtffind4
