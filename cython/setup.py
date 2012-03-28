#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
#
#ext_modules = [Extension("test", ["math.pyx"])]
#setup(
#        name = 'Hello world app',
#        cmdclass = {'build_ext': build_ext},
#        ext_modules = ext_modules
#        )

from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "My hello app",
    ext_modules = cythonize("*.pyx"),
)
