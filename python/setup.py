from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
atmmetaforce_header_dir = '@ATMMETAFORCEPLUGIN_HEADER_DIR@'
atmmetaforce_library_dir = '@ATMMETAFORCEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_atmmetaforce',
                      sources=['ATMMetaForcePluginWrapper.cpp'],
                      libraries=['OpenMM', 'ATMMetaForcePlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), atmmetaforce_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), atmmetaforce_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='atmmetaforce',
      version='0.3',
      py_modules=['atmmetaforce'],
      ext_modules=[extension],
     )
