#!/usr/bin/env python
"""Unit tests for the CNVseq"""
import unittest

import logging
logging.basicConfig(level=logging.ERROR, format="%(message)s")

import warnings
warnings.filterwarnings('ignore', category=ImportWarning)

from cnvsq import reference

class RefTests(unittest.TestCase):
    """Tests for refernce build"""

    pass
