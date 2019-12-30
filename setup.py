import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

__version__ = open("ldetect2/version.py").readline().split(" = ")[1].replace(
    '"', '').strip()
macros = []

setup_requires = ["cython"]
install_requires = [
    "numpy"
]


compile_options = [
    "-Ofast", "-Wall"
]  #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

# print(conda_lib)

extensions = [
    Extension(
        "ldetect2.src.calc_covar", ["ldetect2/src/calc_covar.pyx"],
        language="c"),
    Extension(
        "ldetect2.src.matrix_to_vector", ["ldetect2/src/matrix_to_vector.pyx"],
        language="c"),
    Extension(
        "ldetect2.src.m2v", ["ldetect2/src/m2v.pyx"],
        language="c"),
]


setup(
    name="ldetect2",
    packages=find_packages(),
    ext_modules=cythonize(extensions, annotate=True, language_level='3'),
    # scripts=["bin/epic2", "bin/epic2-df", "bin/epic2-bw"],
    # package_data={
    #     'epic2': [
    #         'effective_sizes/*.txt', 'chromsizes/*chromsizes',
    #         'examples/*.bed.gz'
    #     ],
    #     '': ['*.pyx', '*.pxd', '*.h', '*.c', '*.hpp', "*.bed.gz"]
    # },
    version=__version__,
    description="Fast finding of approx independent linkage disequilibrium blocks in human populations.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/ldetect2",
    keywords=["genetics"],
    license=["BSD"],
    setup_requires=setup_requires,
    install_requires=install_requires,
    include_dirs=["."],
    long_description="."
    )
