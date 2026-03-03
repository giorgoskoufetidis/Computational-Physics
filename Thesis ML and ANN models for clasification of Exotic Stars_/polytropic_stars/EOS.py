import numpy as np 


class EOS:
    """
    In this class we define all the equation of state,
    E = E(P)
    """
    def __init__(self,P):
        self.P = P

    def MDI_1(self):
        return 4.1844 * self.P**0.81449 + 95.00135 * self.P**0.31736
 
    def MDI_2(self):
        return 5.97365 * self.P**0.77374 + 89.24 * self.P**0.30993
    
    def MDI_3(self):
        return 15.55*self.P**0.666 + 76.71 * self.P**0.247
    
    def MDI_4(self):
        return 25.99587 * self.P**0.61209 + 65.62193 * self.P**0.15512
    
    def NLD(self):
        return 119.05736 + 304.80445 * (1 - np.exp(- self.P/48.61465)) + 33722.34448 * (1 - np.exp(- self.P/17499.47411))
    
    def HHJ_1(self):
        return 1.78429 * self.P**0.93761 + 106.93652 * self.P**0.31715
    
    def HHJ_2(self):
        return  1.18961 * self.P**0.96539 + 108.40302 * self.P**0.31264

    def SKa(self):
        return 0.53928 * self.P**1.01394 + 94.31452 * self.P**0.35135
    
    def SkI4(self):
        return 4.75668 * self.P**0.76537 + 105.722 * self.P**0.2745
    
    def HLPS_2(self):
        return 172.858 * (1 - np.exp(- self.P/22.8644)) + 2777.75 * (1 - np.exp(- self.P/1909.97)) + 161.553
    
    def HLPS_3(self):
        return 131.811 * (1 - np.exp(-self.P/4.41577)) + 924.143 * (1 - np.exp(- self.P/523.736)) + 81.5682
    
    def SCVBB(self):
        return 0.371414 * self.P**1.08004 + 109.258 * self.P**0.351019
    
    def WFF_1(self):
        return 0.00127717 * self.P**1.69617 + 135.233 * self.P**0.331471
    
    def WFF_2(self):
        return 0.00244523 * self.P**1.62962 + 122.076 * self.P**0.340401
    
    def PS(self):
        return 9805.95*(1 - np.exp(-0.000193624 * self.P)) + 212.072 * (1 - np.exp(-0.401508 * self.P))+ 1.69483
    
    def W(self):
        return 0.261822 * self.P**1.16851 + 92.4893 * self.P**0.307728
    
    def BGP(self):
        return 0.0112475 * self.P**1.59689 + 102.302 * self.P**0.335526
    
    def BL_1(self):
        return 0.488686 * self.P**1.01457 + 102.26 * self.P**0.355095
    
    def BL_2(self):
        return 1.34241 * self.P**0.910079 + 100.756 * self.P**0.354129
    
    def DH(self):
        return 39.5021 * self.P**0.541485 + 96.0528 * self.P**0.00401285

    def APR_1(self):
        return 0.000719964 * self.P**1.85898 + 108.975 * self.P**0.340074


    def models(self):
        return [self.MDI_1(),self.MDI_2(),self.MDI_3(),self.MDI_4(),self.NLD(),self.HHJ_1(),self.HHJ_2(),self.SKa(),
                self.SkI4(),self.HLPS_2(),self.HLPS_3(),self.SCVBB(),self.WFF_1(),self.WFF_2(),self.PS(),self.W(),self.BGP(),self.BL_1(),
                self.BL_2(),self.DH(),self.APR_1()]

