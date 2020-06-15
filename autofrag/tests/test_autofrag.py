"""
Unit and regression test for the autofrag package.
"""

# Import package, test suite, and other packages as needed
import autofrag
import pytest
import sys

def test_autofrag_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "autofrag" in sys.modules
