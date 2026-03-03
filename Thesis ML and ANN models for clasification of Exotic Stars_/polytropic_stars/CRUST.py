import numpy as np 


class CRUST:
    
    def crust_equation1(self, P):
        return 103.17338*(1 - np.exp(-P/0.38527)) + 7.34979* (1 - np.exp(-P/0.01211))+ 0.00873

    def crust_equation2(self, P):
        return 0.00015 + 0.00203 * (1 - np.exp(-P * 344827.5)) + 0.10851 * (1 - np.exp(-P * 7692.3076))

    def crust_equation3(self, P):
        return 0.0000051 * (1 - np.exp(- P * 0.2373 * 1e10) ) + 0.00014 * (1 - np.exp(- P * 0.4020 * 1e8))

    def crust_equation4(self, P):
        c0 = 31.93753
        c1 = 10.82611 * np.log(P) / np.log(10)
        c2 = 1.29312 * (np.log(P) / np.log(10))**2
        c3 = 0.08014 * (np.log(P) / np.log(10))**3
        c4 = 0.00242 * (np.log(P) / np.log(10))**4
        c5 = 0.000028 * (np.log(P) / np.log(10))**5
        return 10**(c0 + c1 + c2 +c3 + c4 + c5)
    def equation(self, P):
        if 9.34375 * 1e-5 <= P <= 0.184:
            return self.crust_equation1(P)

        if 4.1725 * 1e-8  <= P < 9.34375 * 1e-5:
            return self.crust_equation2(P)

        if 1.44875 * 1e-11 <= P < 4.1725 * 1e-8:
            return self.crust_equation3(P)

        if P < 1.44875 * 1e-11:
            return self.crust_equation4(P)
