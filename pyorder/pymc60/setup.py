#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pymc60', parent_package, top_path)


    config.add_extension(
        name='mc60module',
        sources=['src/mc60.pyf', 'src/mc60ad.f'],
        libraries=[],
        library_dirs=[],
        include_dirs=['src'],
        extra_link_args=[])

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        #package_dir={'pymc60':'Lib'},
        #packages = ['pymc60'],
        **configuration(top_path='').todict())
