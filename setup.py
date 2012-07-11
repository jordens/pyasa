#!/usr/bin/python
# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Lester Ingber's Adaptive Stimulated
#   Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@phys.ethz.ch>

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy
import os

setup(
    name="pyasa",
    description="python bindings for Lester Ingber's Adaptive Stimulated Annealing",
    long_description= """http://www.ingber.com/#ASA""",
    version="0.0+dev",
    author="Robert Jordens",
    author_email="jordens@gmail.com",
    url="http://launchpad.net/pyasa",
    license="modified BSD",
    install_requires=["numpy"],
    #packages=["flycapture2"],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("asa",
        define_macros = [('OPTIONAL_DATA_PTR', '1'),
                         ('USER_TYPE', 'void *'),
                         ('ASA_TEMPLATE', '0'),
                         ('ASA_LIB', '1')],
        extra_compile_args=["-g"],
        sources = ["src/asa.pyx", "ASA/asa.c", "src/asa_rand.c"],
        #, "ASA/asa_usr.c"],
        include_dirs = ["ASA", numpy.get_include(),],
        ),]
