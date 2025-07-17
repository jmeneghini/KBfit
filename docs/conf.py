import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath('..'))

# Project information
project = 'KBfit'
author = 'KBfit Development Team'
copyright = f'{datetime.now().year}, {author}'
version = '1.0'
release = '1.0.0'

# Extensions
extensions = [
    'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
]

# Try to import myst_parser
try:
    import myst_parser
    extensions.append('myst_parser')
except ImportError:
    print("Warning: myst_parser not found, Markdown support disabled")
    # Remove .md files from source_suffix if myst_parser is not available
    source_suffix = {'.rst': None}

# Breathe configuration for Doxygen integration
breathe_projects = {'KBfit': os.path.join(os.path.dirname(__file__), '_build', 'doxygen', 'xml')}
breathe_default_project = 'KBfit'
breathe_default_members = ('members', 'undoc-members')

# MyST parser configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# HTML theme configuration
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# HTML output configuration
html_title = f'{project} Documentation'
html_short_title = project
html_show_sourcelink = True
html_show_sphinx = True
html_show_copyright = True
html_last_updated_fmt = '%b %d, %Y'
html_use_index = True
html_split_index = False
html_copy_source = True
html_show_sourcelink = True
html_sourcelink_suffix = ''
html_use_opensearch = ''
html_file_suffix = None
html_link_suffix = None

# LaTeX output configuration
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'fncychap': '',
    'maketitle': '',
    'printindex': '',
    'sphinxsetup': '',
}

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# TODO configuration
todo_include_todos = True
todo_emit_warnings = True

# Math configuration
mathjax_config = {
    'tex2jax': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
        'processEscapes': True,
        'processEnvironments': True,
    },
    'displayAlign': 'center',
    'displayIndent': '0em',
}

# Source file configuration
if 'myst_parser' in extensions:
    source_suffix = {
        '.rst': None,
        '.md': 'myst_parser',
    }
else:
    source_suffix = {'.rst': None}

# Master document
master_doc = 'index'
