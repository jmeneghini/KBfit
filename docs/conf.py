import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'KBfit'
author = 'KBfit developers'

extensions = ['breathe', 'myst_parser']

breathe_projects = {'KBfit': os.path.join(os.path.dirname(__file__), '_build', 'doxygen', 'xml')}
breathe_default_project = 'KBfit'

html_theme = 'sphinx_rtd_theme'
