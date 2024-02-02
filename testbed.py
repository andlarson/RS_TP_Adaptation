"""
This file acts as a testbed for individual components of the library.
"""


class MyC:
    
    def __init__(self, a, b):
        
        self.a = a
        self.b = b

    def __setattr__(self, name, val):
        
        if name == "a":
            print("The a var is being set....But I'm not going to actually set it.")

        if name == "b":
            print("The b var is being set...And I'm actually going to set it.")
            super().__setattr__(name, val)


if __name__ == "__main__":
    
    c = MyC(10, 20)

    print(c.b)

