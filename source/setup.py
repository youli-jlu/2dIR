#from setuptools import setup
from setuptools import Extension, setup
#from distutils.extension import Extension

from Cython.Build import cythonize

libraries = [ "fftw3","m"]
extra_compile_args = ['-O3', '-std=c++14']
extra_link_args = ['-O3', '-std=c++14']

ext = Extension(
        "spectra2d",
        ["spectra2d.pyx",],
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        )

#setup(ext_modules=cythonize("spectra2d.pyx"))
setup(ext_modules=cythonize([ext,]))
