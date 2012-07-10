#!/usr/bin/python
# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Lester Ingber's Adaptive Stimulated
#   Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@phys.ethz.ch>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    license="hybrid BSD",
    install_requires=["numpy"],
    #packages=["asa"],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("asa",
        define_macros = [('OPTIONAL_DATA_PTR', '1'),
                         ('USER_TYPE', 'void *'),
                         ('ASA_LIB', '1')],
        extra_compile_args=["-g"],
        sources = ["src/asa.pyx", "ASA/asa.c", "src/asa_rand.c"],
        include_dirs = ["ASA", numpy.get_include(),],
        ),]
)
