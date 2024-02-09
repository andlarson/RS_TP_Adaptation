"""
This file acts as a testbed for individual components of the library.
"""


class LL_entry:

    def __init__(self, prev: "LL_entry", val: int):
        self.prev = prev
        self.val = val



class LL:
    
    def __init__(self):
        self.l = []

    def add(self, val):
        if len(self.l) == 0:
            entry = LL_entry(None, val)
        else:
            prev = self.l[-1]
            entry = LL_entry(prev, val)
        return entry

    def remove(self, entry):
        

    


if __name__ == "__main__":
    
