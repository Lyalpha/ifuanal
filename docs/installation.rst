Installation
============

Requirements
------------

ifuanal has the following python dependencies:

* python 2.7.11
* numpy 1.10
* scipy 0.15
* matplotlib 1.3
* astropy 1.1.1
* astroquery 0.2.3 (optional: allows automated retrieval of Galactic
  extinction) 
* dill 0.2.5
* pathos 0.2

A recommended way to set this up is to install `astroconda
<http://astroconda.readthedocs.io/en/latest/installation.html>`_ and
install extra dependancies from within your astroconda environment
with::

     $ pip install astroquery
     $ pip install dill
     $ pip install pathos

.. _starlight-install:

STARLIGHT installation
----------------------

The STARLIGHT executable suitable for your system should be downloaded from
`here <http://astro.ufsc.br/starlight/node/3>`_. This should be copied (or linked) in to the `starlight/` subdirectory of the ifuanal project directory with the exact file name ``StarlightChains_v04.exe`` so that the code can find it.

ifuanal comes bundled with a set of model bases from Bruzual & Charlot (2003), although others can be added manually (see the :ref:`Tutorial <cont-fitting>`).

Setup
-----
No installation is necessary, just be sure to add the project directory
directory to your `PATH` if you wish to import the module from
anywhere. Otherwise, from within the project directory::

    >>> import ifuanal

See :doc:`the tutorial <tutorial>` for example workflows.
