from setuptools import setup
from glob import glob

import hts.htsffi
import cffi

# from mpld3
def get_version(fname):
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(fname) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")


#cffi.verifier.cleanup_tmpdir()
ext = hts.htsffi.ffi.verifier.get_extension()

setup(name='hts',
      version=get_version("hts/__init__.py"),
      packages=['hts'],
      package_data={'hts': ["hts/hts_concat.h", "hts/hts_extra.h", "hts/hts_extra.c"]},
      ext_package='htsffi',
      include_package_data=True,
      ext_modules=[ext],
      zip_safe=False,
      install_requires=['cffi', 'setuptools>=0.6c11']
      )

