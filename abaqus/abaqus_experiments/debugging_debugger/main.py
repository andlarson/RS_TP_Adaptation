"""
This file is used to test the debugger.
"""

import sys
sys.path.insert(0, "/home/andlars/Desktop/RS_TP_Adaptation/experiments/experiments/debugging_debugger/")

import test


if __name__ == "__main__":
    
    a = 5
    b = 10
    print(str(a+b))
    
    test.test_function()
