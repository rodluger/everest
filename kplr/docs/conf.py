# -*- coding: utf-8 -*-

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Sphinx setup.
extensions = ["sphinx.ext.autodoc", ]  # 'sphinx.ext.mathjax']
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"

# Project specs.
project = u"kplr"
copyright = u"2013, Dan Foreman-Mackey"
version = "0.1.0"
release = "0.1.0"

# General options.
exclude_patterns = ['_build']
pygments_style = 'sphinx'

# Theme setup.
html_theme = "kplr"
html_theme_path = ["_themes"]
# html_theme_options = {}
html_static_path = ['_static']
html_use_smartypants = True
html_sidebars = {
    "**": ["relations.html"]
}
# html_additional_pages = {}
htmlhelp_basename = "kplrdoc"
