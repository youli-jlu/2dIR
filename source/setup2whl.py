from setuptools import Extension, setup
#from distutils.extension import Extension

from Cython.Build import cythonize
#from Cython.Distutils import build_ext
from Cython.Build import build_ext

libraries = [ "fftw3","m"]
extra_compile_args = ['-O3', '-std=c++14']
extra_link_args = ['-O3', '-std=c++14']
print("get the flags")

ext = Extension(
        #"spectra1d",
        #["spectra1d.pyx",],
        "spectra2d",
        ["spectra2d.pyx",],
        libraries=libraries,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        )
print("create the ext")

setup(name="spectra2d",
        version="3.0",
        ext_modules=cythonize(
            [ext,],
            build_dir="./",
            compiler_directives=dict(always_allow_keywords=True)
            ),
        cmdclass=dict( build_ext=build_ext),packages=[]
        )

#setup(
#    name="MyModule"  # 模块名称 import MyModule,
#    ext_modules=cythonize(
#        [
#           Extension("pkg1.*", ["root/pkg1/*.py"]),
#           Extension("pkg2.*", ["root/pkg2/*.py"]),
#           Extension("1.*", ["root/*.py"])
#        ],
#        build_dir="build",
#        compiler_directives=dict(
#        always_allow_keywords=True
#        )),
#    cmdclass=dict(
#        build_ext=build_ext
#    ),
#    packages=["pkg1", "pkg2"]  # packages=[]时打包后的wheel文件中不含源码(.py)
#)
