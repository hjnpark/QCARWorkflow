"""
setup.py: Install QCArchive Workflow script.  
"""
import versioneer
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='qcarw',
    description='An automated workflow that can refine MD simulation trajectories.',
    url='https://github.com/hjnpark/QCARWorkflow',
    author='Heejune Park',
    packages=find_packages(),
    package_data={'': ['*.ini']},
    include_package_data=True,
    long_description=long_description,
    long_description_content_type="text/markdown",
    scripts=['bin/Nebterpolate.py'],
    entry_points={'console_scripts': [
        'qcarw-refine = qcarw.refine:main',
        'qcarw-dsoptimize = qcarw.dsopt:main',
        'qcarw-smooth = qcarw.smooth:main',
        'qcarw-optimize = qcarw.opt:main',
        'qcarw-neb = qcarw.neb:main',
    ]},
    install_requires=[
        'numpy>=1.11',
        'networkx',
        'scipy>=1.6',
        'matplotlib',
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
    ],
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
    zip_safe=True,
)

