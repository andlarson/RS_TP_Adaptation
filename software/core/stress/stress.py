import util.geom as geom



# Compute how near the stress field is to equilibrium.
def check_equil():
    pass



class StressTensor:
    
    def __init__(self, s11, s22, s33, s13, s12, s23):
    # type: (float, float, float, float, float, float) -> None

        self.s11 = s11
        self.s22 = s22
        self.s33 = s33
        self.s13 = s13
        self.s23 = s23
        self.s12 = s12



class StressProfile:

    def __init__(self, stress_profile_list):
    # type: (list[tuple[geom.Point3D, StressTensor]]) -> None 

        # TODO: Check how close to equilibrium the stress field is. Warn the
        #   user if the stress field is not in equilibrium.

        self.stress_profile_list = stress_profile_list
