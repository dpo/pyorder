#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    import numpy
    import ConfigParser
    import os
    from numpy.distutils.misc_util import Configuration

    # Read relevant PyOrder-specific configuration options.
    pyorder_config = ConfigParser.SafeConfigParser()
    pyorder_config.read(os.path.join(top_path, 'site.cfg'))
    hsl_dir = pyorder_config.get('HSL', 'hsl_dir')

    config = Configuration('pymc21', parent_package, top_path)

    config.add_extension(
        name='mc21module',
        sources=['src/mc21.pyf', os.path.join(hsl_dir,'mc21d','mc21d.f')],
        libraries=[],
        library_dirs=[],
        include_dirs=['src'],
        extra_link_args=[])

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
