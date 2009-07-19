#!/usr/bin/env python
"""
PyOrder: Sparse Matrix Ordering Methods in Python

PyOrder is a library of ordering methods to reduce the profile, wavefront or
semi-bandwidth or symmetric and non-symmetric sparse matrices.
D. Orban <dominique.orban@gerad.ca>
"""

DOCLINES = __doc__.split("\n")

import os
import sys

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: LGPL
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('pyorder')

    # Set config.version
    config.get_version(os.path.join('pyorder','version.py'))

    return config

def setup_package():

    from numpy.distutils.core import setup, Extension
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'pyorder')) # to retrieve version

    try:
        #extra_link_args=[]
        #pymc60 = Extension(
        #             name='mc60module',
        #             sources=['Src/mc60.pyf', 'Src/mc60ad.f'],
        #             libraries=[],
        #             library_dirs=[],
        #             include_dirs=['Src'],
        #             extra_link_args=extra_link_args)

        setup(
            name = 'pyorder',
            author = "Dominique Orban",
            author_email = "dominique.orban@gerad.ca",
            maintainer = "PyOrder Developers",
            maintainer_email = "dominique.orban@gerad.ca",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "",
            download_url = "",
            license = 'LGPL',
            classifiers=filter(None, CLASSIFIERS.split('\n')),
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            configuration=configuration,
            #ext_modules = [pymc60],
            #package_dir={'pyorder':'Lib'},
            #packages = ['pyorder'],
            )
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
