# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
# with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='pyglib',
    version='1.0.0',
    description='Python Modules for CyGutz package',
    # The project's main homepage.
    # url='https://github.com/pypa/sampleproject',
    # Author details
    # author='Yongxin Yao, Nicola Lanata',
    # author_email='cygutz@gmail.com',
    license='OpenBSD',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: OpenBSD License',

        # Specify the Python versions you support here.
        'Programming Language :: Python :: 2.7',
    ],

    #packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    packages=["pyglib"],

    install_requires=['numpy', 'scipy'],

)
