
from __future__ import print_function

import os

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

this_directory = os.path.dirname(__file__)
source_directory = os.path.join(this_directory, 'chull')
exec(open(os.path.join(source_directory, 'version.py')).read())  # Load in the variable __version__.

dependencies = ['cypari']

setuptools.setup(
    name='chull',  
    version=__version__,
    author="William Worden",
    author_email="wtworden@gmail.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wtworden/CHull",
    packages=['chull'],
    package_data={
        },
    install_requires=dependencies,
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
 )


## sage -python setup.py sdist bdist_wheel --universal
## twine upload dist/*