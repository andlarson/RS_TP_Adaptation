"""
Utilities useful for debugging.
"""

import sys


def dp(message: str) -> None:
    """Prints a message to the standard error stream."""

    sys.__stderr__.write(message + "\n")
