"""
Utilities useful for debugging.
"""

import sys
import traceback
import math



def dp(message: str) -> None:
    """Prints a message to the standard error stream."""

    sys.__stderr__.write(message + "\n")



def dump_exception() -> None:
    """Dumps information related to the most recent exception which occurred."""

    # Dump the exception being handled. 
    exception_info = sys.exc_info()
    
    dump_banner("EXCEPTION OCCURRED")
    dp("")
    dp("TYPE: " + str(exception_info[0]))
    dp("")
    dp("CAUSE: " + str(exception_info[1]))
    dp("") 
    dp("TRACEBACK: ")
    tb = exception_info[2]
    traceback.print_tb(tb, file=sys.__stderr__)
    dp("")
    dump_banner_end()



BANNER_LEN = 80
STAR_CNT = 4
def dump_banner(msg: str) -> None:
    """Prints a message embedded in a formatted banner."""
    
    stars = STAR_CNT * "*" 

    dp("")
    dp(BANNER_LEN * "*")

    room_for_whitespace = BANNER_LEN - len(msg) - 2 * STAR_CNT
    whitespace_per_side = room_for_whitespace / 2
    if room_for_whitespace % 2 == 1:
        # Odd number of characters available for whitespace, so uneven padding.
        whitespace_buffer1 = math.floor(whitespace_per_side) * " "
        whitespace_buffer2 = math.ceil(whitespace_per_side) * " "
        dp(stars + whitespace_buffer1 + msg + whitespace_buffer2 + stars)
    else:
        # Even number of characters available for whitespace, so even padding.
        whitespace_buffer = int(whitespace_per_side) * " "
        dp(stars + whitespace_buffer + msg + whitespace_buffer + stars)





def dump_banner_end() -> None:
    """Prints the end of a banner."""
    
    spaces = BANNER_LEN - 2 * STAR_CNT
    whitespace = spaces * " "
    stars = STAR_CNT * "*"

    dp(stars + whitespace + stars)
    dp(BANNER_LEN * "*")
    dp("")
