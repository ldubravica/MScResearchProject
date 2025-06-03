from setuptools import setup
from Cython.Build import cythonize

setup(
    name='lz76',
    ext_modules=cythonize("lz76.pyx"),
)