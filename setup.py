


from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys

#gmpgcd = Extension("gmpgcd", include_dirs=['/opt/local/include','usr/local/include'], libraries=['gmp'], sources=['gmpgcdmodule.c'])

static_libraries = ['gmp']
static_lib_dir = '/usr/local/lib'

if sys.platform == 'win32':
    libraries.extend(static_libraries)
    library_dirs.append(static_lib_dir)
    extra_objects = []
else: # POSIX
    extra_objects = ['{}/lib{}.a'.format(static_lib_dir, l) for l in static_libraries]

gmpgcd = Extension("gmpgcd",include_dirs=['usr/local/include'],extra_objects=extra_objects,sources=['gmpgcdmodule.c'],extra_link_args=['--disable-assembly'])
                

setup(
    name="gmpgcd",
    version='1.0',
    description='GMP computation of gcd',
    ext_modules=[gmpgcd]
)