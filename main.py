# Implementation of FFT
import math
import numpy as np
from cmath import pi

#p(x) = x + 2x^2 + 3x^3
#g(x) = 3 + 5x + 4x^4

# Degree 7 : need 8 pts
# evaluate polynomials at 8 roots of unity

# Custom Complex Class
class Complex:
    def __init__(self,a,b=0):
        self.a = a
        self.b = b
    def __repr__(self):
        char = "0" if abs(round(self.b)) == 0 else str(abs(round(self.b,2)))
        im = "i" if char == "1" else char+"i"
        re = "0" if round(self.a) == 0 else str(round(self.a,2))
        if(self.b >= 0):
            if re == "0" and im == "0i":
                s = "0"
            elif re == "0":
                s = im
            elif im == "0i":
                s = re
            else:
              s = f"({re} + {im})"  
        else:
            if re == "0" and im == "0i":
                s = "0"
            elif re == "0":
                s = im
            elif im == "0i":
                s = re
            else:
              s = f"({re} - {im})"  
        return s
    def __add__(self,other):
        return Complex(self.a+other.a,self.b+other.b)
    def __mul__(self,other):
        #a + ib * c + id = ac +i(ad + bc) - bd
        return Complex(self.a * other.a -self.b*other.b,self.a*other.b + self.b*other.a)
    def __sub__(self,other):
        return Complex(self.a-other.a,self.b-other.b)
    def __abs__(self):
       mod = self.a **2 + self.b** 2
       return math.sqrt(mod)
    def conjugate(self):
        return Complex(self.a,-self.b)
    def __truediv__(self,other):
        if isinstance(other,float):
            return Complex(self.a/other,self.b/other)
        else:
            denom = float(other.__abs__() ** 2)
            conjugate = other.conjugate()
            num = self * conjugate
            # print(num)
            return num/denom
    def __pow__(self,power):
        if power == 0 :
            return Complex(1)
        return self * self.__pow__(power - 1)
    
    @classmethod
    def nthroots(self,n):
        c = 2*pi/n
        return Complex(math.cos(c),math.sin(c))
        # return Complex(self * other.conjugate()/denom)
        
#Run time: O(nlogn)
def FFT(poly):
    # P = [p0,p1,....pn-1] # degree n -1
    # From coeff to value repr(sample space that is graphical)
    n = len(poly) # Degree of polynomial power of 2
    if n == 1: # means degree 0 that is constant
        return poly
    # omega = np.exp(2j*np.pi/n)
    omega = Complex.nthroots(n)
    p_e = poly[::2]
    p_o = poly[1::2]
    y_e,y_o = FFT(p_e),FFT(p_o)

    y = [Complex(0,0)] * n
    for j in range(int(n/2)):
        y[j] = y_e[j] + pow(omega,j)*y_o[j]
        y[j+n//2] = y_e[j] - pow(omega,j)*y_o[j]
    return y

#Run time: O(nlogn)
def IFFT(poly):
    def IFFT_pre(poly):
    # From value rep(Sample space that is graph) to coeff
        n = len(poly) # Degree of polynomial power of 2
        if n == 1: # means degree 0 that is constant
            return poly
        # omega = np.exp(2j*np.pi/n)
        omega = Complex.nthroots(-n)
        p_e = poly[::2]
        p_o = poly[1::2]
        y_e,y_o = IFFT_pre(p_e),IFFT_pre(p_o)
        y = [Complex(0,0)] * n
        for j in range(int(n/2)):
            y[j] = y_e[j] + pow(omega,j)*y_o[j]
            y[j+n//2] = y_e[j] - pow(omega,j)*y_o[j]
        return y
    l = IFFT_pre(poly)
    return [item / Complex(len(poly)) for item in l]

def multiply(poly1,poly2):
    p1 = FFT(poly1)
    p2 = FFT(poly2)
    if len(p1) != len(p2):
        if len(p1) < len(p2):
            for i in range(len(p2) - len(p1)):
                p1.append(Complex(0))
        else:
            for i in range(len(p1) - len(p2)):
                p2.append(Complex(0))
    print(p1)
    print(p2)
    l = []
    for i in range(len(p1)):
        l.append(p1[i] * p2[i])
    return IFFT(l)
p1 = [Complex(0),Complex(1),Complex(2),Complex(0),Complex(0),Complex(0)] # x + 2x^2
p2 = [Complex(1),Complex(0),Complex(1),Complex(3),Complex(0),Complex(0)] # 1 + x^2 + 3x^3

y1 = FFT(p1)
y2 = FFT(p2)
print(y1)
print(y2)

print(multiply(p1,p2))