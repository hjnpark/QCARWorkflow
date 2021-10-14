"""
setup.py: Install QCArchive Workflow script.  
"""

VERSION = 1.1
__author__ = "Heejune Park"
__version__ = "%.1f"%VERSION

import os, sys
from distutils.core import setup#, Extension
import numpy
import glob

requirements = ['numpy', 'networkx', 'qcelemental', 'qcfractal', 'qcportal', 'qcengine', 'geometric', 'psi4']
    
def buildKeywordDictionary():
    setupKeywords = {}
    setupKeywords["name"]              = "QCARWorkflow"
    setupKeywords["version"]           = "%.1f" %VERSION
    setupKeywords["author"]            = __author__
    setupKeywords["author_email"]      = "heepark@ucdavis.edu"
    setupKeywords["packages"]          = ["QCARWorkflow", "nebterpolator", "nebterpolator.io", "nebterpolator.core"]
    setupKeywords["package_dir"]       = {"QCARWorkflow": "src"}
    setupKeywords["scripts"]           = glob.glob("bin/*.py")# + glob.glob("bin/*.sh") + glob.glob("bin/*.exe") + glob.glob("bin/*.vmd")
    setupKeywords["platforms"]         = ["Linux"]
    setupKeywords["description"]       = " An automated workflow that can refine reaction pathways from MD simulation trajectories."
    outputString=""
    firstTab     = 40
    secondTab    = 60
    for key in sorted( setupKeywords.keys() ):
         value         = setupKeywords[key]
         outputString += key.rjust(firstTab) + str( value ).rjust(secondTab) + "\n"
    print("%s" % outputString)
    return setupKeywords
    

def main():
    setup_keywords = buildKeywordDictionary()
    setup(**setup_keywords)
    for requirement in requirements:
      try:
          exec('import %s' % requirement)
      except ImportError as e:
          print('\nWarning:%s' % e, file=sys.stderr)
          print('Warning: Some package functionality may not work', file=sys.stderr)

if __name__ == '__main__':
    main()

