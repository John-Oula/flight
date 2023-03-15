

class UAV(object):
    """
    Initialize UAV configuration and flight conditions

    """
    def __init__(self,b,S,h,V,c,theta_zero,m,I_xx,I_yy,I_zz,I_xz):
        self.wing_span = b
        self.wing_area = S
        self.altitude = h
        self.velocity = V
        self.mean_ac = c
        self.theta_0 = theta_zero
        self.mass = m
        self.I_xx = I_xx
        self.I_yy = I_yy
        self.I_zz = I_zz
        self.I_xz = I_xz


