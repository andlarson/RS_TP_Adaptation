"""
This file acts as a testbed for individual components of the library.
"""

import copy


class MyC:
    
    def __init__(self, a, b, c):
        
        self.a = a
        self.b = b
        self.dict = {"A": MyOtherC(c, 20), "B": MyOtherC(c, 700)}


class MyOtherC:
    
    def __init__(self, e, f):
        self.e = e
        self.f = f 


if __name__ == "__main__":
    
    c = MyC(10, 20, 300)
    c.dict["D"] = MyOtherC(987, 65)
    print(c.dict)

    d = copy.deepcopy(c)
    print(d.dict)
    print(d.dict["A"].e)
    print(d.dict["B"].e)
