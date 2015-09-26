# Michael Hentrich is altering this script so as to create pedagogical
# version which will display each step of the process.
# Date: September 2015
#
# xed(a,b) for displaying steps of the process
# xlgcd(a,b) - like xed but includes detailed summary
# isolve(a,b,d) - solves Diophntine equations using euclid
# iterative_egcd(a, b)
# modinv(a, m)
# icrt(a, n) - Chinese Remainder Theorem
# mul_inv(a, b) - multiplicative modular inverse
#
# The MotherLoad for Number Theory: see Brian Gladman
# 173.254.28.24/~brgladma/number_theory.py
# crt, gcd, xgcd, lcm, cong, inv
# 
# Kirby Urner for a cooler version of the Extended Euclidean Algorithm, with
# related inverse and Chinese Remainder Theorem functions.
# euclid(a,b)
# ext_euclid(a,b)
#
# Functions at the the end are curosities where I am testing out
# using product, add, subtract from numpy to use with reduce()
# eea_map, xcrt, inverse
#
# Author: Sam Erickson for test(a,b) and getCoeffic(a,b)
# Date: 7/14/2013
#
# Program Description: This program gives the integer coefficients x,y to the
# equation ax+by=gcd(a,b) given by the extended Euclidean Algorithm. Note that
# this program does not use the typical extended Euclidean Algorithm to find
# ax+by=gcd(a,b).
"""
General algorithm for program:
1. Get integers a and b from user.

2. Call test(a,b) and assign to get one integer coefficient for the equation
(1/gcd(a,b))*[a*x+b*y=gcd(a,b)].


3. Call getCoffic(min2,aNew,bNew,GCD) to return a string for of the equation
a*x+b*y=gcd(a,b) for printing.
"""

def test(a,b):
    """
    Input: a and b must be integers.

    Output: If solutionsDict isn't empty, then the min positive integer along
    with a/gcd(a,b), b/gcd(a,b), and gcd(a,b) will be returned. If solutionsDict
    does not contain positive integers the max negative integer will be returned
    instead of the least positive integer.
    Otherwise if solutionsDict is empty, None,None,None,None will be returned.
    """
    from fractions import gcd
    # Test which of a and b is larger.
    if b>a:
        a,b=b,a
    
    solutionsDict={}

    # Iterate range(1,aNew+1) to find solutions to [1-bNew*x]%aNew==0
    for x in range(1,abs(a//gcd(a,b))):
        if (1-(b//gcd(a,b))*x)%(a//gcd(a,b))==0:
            solutionsDict[x]=x
    # Determine if solutionsDict is empty. 
    if len(solutionsDict)!=0:
        # Make sure solutionsDict contains positive integers
        ct=0
        for x in solutionsDict:
            if x>0:
                ct+=1
                break
        
        # If solutionsDict contains positive integers get least positive integer
        if ct>0:
            while min(solutionsDict)<=0:
                del d1[min(solutionsDict)]
    
            return min(solutionsDict),a//gcd(a,b),b//gcd(a,b),gcd(a,b)

        else:
            return max(solutionsDict),a//gcd(a,b),b//gcd(a,b),gcd(a,b)
            
    else:
        return None,None,None,None

def getCoeffic(Min,aPrime,bPrime,GCD):
    """
    Input: Min,aPrime,bPrime,GCD is the 4-tuple obtained from test(a,b), in the
    same order.

    Output: If the final expression ax+by really does equal gcd(a,b),and if
    Min isn't None, the string "ax+by=gcd(a,b)" will be returned. Otherwise
    if Min=None or output string "ax+by" is not equal to gcd(a,b), error
    messages will be printed accordingly.
    """

    # Create string for ax+by=gcd(a,b)
    if Min!=None:
        cA=int((1-Min*bPrime)/aPrime)
        cB=Min
        s1=str(cA)+"*"+str(aPrime*GCD)+" + "+str(cB)+"*"+str(bPrime*GCD)
        s2="="+str(GCD)
        
        # Test that string expression is correct before returning it 
        if eval(s1)==GCD:
            return s1+s2
        # Print Error Message
        else:
            print("The program must have a bug since in the final expression ax+by was not equal to gcd(a,b).")
    # Print Error Message
    else:
        print("Either program bug caused no results or algorithm failed.")

def euclid(a, b):
    """Euclid's algorithm for GCD

    Given input a, b the function returns d such that gcd(a,b) = d"""

    if a < b:
        a,b = b,a
    else:
        pass

    while b != 0:
        a, b = b, a % b
    return a

def ext_euclid(a, b):
    """Extended Euclid's algorithm for GCD

    Given input a, b the function returns d such that gcd(a,b) = d
    and x, y such that ax + by = d, as well as u,v such that au=bv"""

    if a < b:
        a,b = b,a
    else:
        pass
    u, v, x, y = 0, 1, 1, 0
    while b != 0:
        a, b, x, y, u, v = b, a%b, u, v, x - (a // b) * u, y - ( a // b) * v
    return a, x, y, u, v

def xed(a, b):
    """Displaying the steps of Extended Euclid's algorithm for GCD

    Given input a, b the function returns d such that gcd(a,b) = d
    and x, y such that ax + by = d, as well as u,v such that
    xv + yu = d"""

    d = euclid(a, b)
    step = 0
    
    if a < b:
        a,b = b,a
    else:
        pass
    u, v, x, y = 0, 1, 1, 0

    A, B = a, b

    print "Displaying the steps of Extended Euclid's algorithm for GCD\n"
    print "The intitial values are:"
    print "a = ", a, "\nb = ", b
    print "x = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
    print "----------------------------------\n"
    print "We write the equation in this form: a = q*b + r"
    print "\n", a, "=", a//b,"*", b, "+", a%b
    print "\nwhere q is the integer part of a/b, and r is a%b\n"
    print "----------------------------------\n"
    print "During the recursive loop, values will change.\n"
    print "a <---| b"
    print "b <---| a%b"
    print "x <---| u"
    print "y <---| v"
    print "u <---| x - (a/b)*u"
    print "v <---| y - (a/b)*v"
    print "----------------------------------\n"

    while b != 0:
        a, b, x, y, u, v = b, a%b, u, v, x - (a // b) * u, y - ( a // b) * v
        step += 1
        print "\nSTEP#: ", step
        print "a = ", a, "\nb = ", b
        print "x = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
        print "\nNow,", b, "=", u,"*", A,"+", v,"*", B
        if (b != 0):
            print "Meanwhile,", a, "=", a//b,"*", b, "+", a%b
   
    print "\nx = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
    print "\n", A,"*", x," + ", B,"*", y," = ", (A*x + B*y)
    print "\nDoes", A,"*", x," + ", B,"*", y," = ", d," ?"
    print (A*x + B*y) == d
    print "--------------------\n"

    if ((A*x + B*y) == d):
        print "a, x, y, u, v: "
        return a, x, y, u, v
    else:
        print "There's been some kind of error."

def xlgcd(a, b):
    """Displays the steps of Extended Euclid's algorithm for GCD
    This version displays a summary before returning value.
    It stores these values in lists throughout the loop.


    Given input a, b the function returns d such that gcd(a,b) = d
    and x, y such that ax + by = d, as well as u,v such that
    xv + yu = d"""

    d = euclid(a, b)
    step = 0
    
    if a < b:
        a,b = b,a
    else:
        pass
    u, v, x, y = 0, 1, 1, 0

    La=[]
    Lb=[]
    Lx=[]
    Ly=[]
    Lu=[]
    Lv=[]

    A, B = a, b

    print "Displaying the steps of Extended Euclid's algorithm for GCD\n"
    print "The intitial values are:"
    print "a = ", a, "\nb = ", b
    print "x = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
    print "----------------------------------\n"
    print "We write the equation in this form: a = q*b + r"
    print "\n", a, "=", a//b,"*", b, "+", a%b
    print "\nwhere q is the integer part of a/b, and r is a%b\n"
    print "----------------------------------\n"
    print "During the recursive loop, values will change.\n"
    print "a <---| b"
    print "b <---| a%b"
    print "x <---| u"
    print "y <---| v"
    print "u <---| x - (a/b)*u"
    print "v <---| y - (a/b)*v"
    print "----------------------------------\n"

    La.append(a)
    Lb.append(b)
    Lx.append(x)
    Ly.append(y)
    Lu.append(u)
    Lv.append(v)
    
    while b != 0:
        a, b, x, y, u, v = b, a%b, u, v, x - (a // b) * u, y - ( a // b) * v
#        a = b
#        b = a%b
#        x = u
#        y = v
#        u = x - (a//b) * u
#        v = y - (a//b) * v
        step += 1
        La.append(a)
        Lb.append(b)
        Lx.append(x)
        Ly.append(y)
        Lu.append(u)
        Lv.append(v)
        print "\nSTEP#: ", step
        print "a = ", a, "\nb = ", b
        print "x = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
        print "\nNow,", b, "=", u,"*", A,"+", v,"*", B
        if (b != 0):
            print "Meanwhile,", a, "=", a//b,"*", b, "+", a%b
   
    print "\nx = ", x, "\ny = ", y, "\nu = ", u, "\nv = ", v
    print "\n", A,"*", x," + ", B,"*", y," = ", (A*x + B*y)
    print "\nDoes", A,"*", x," + ", B,"*", y," = ", d," ?"
    print (A*x + B*y) == d
    print "--------------------\n"

    if ((A*x + B*y) == d):
        print "\nSUMMARY of:"
        print "a = ((a/b) * b) + a%b"
        print "------------------------------------"
        for i in range(step):
            print La[i], "= (", La[i]//Lb[i],"*", Lb[i], ") +", La[i]%Lb[i]
        print "\nSUMMARY of:"
        print "r = (u * a0)  +  (v * b0)"
        print "------------------------------------"
        for i in range(step):
            print Lb[i], "= (", Lu[i],"*", A,") + (", Lv[i],"*", B, ")"
        print "\na, x, y, u, v: "
        return a, x, y, u, v
    else:
        print "There's been some kind of error."
    raw_input("Press enter")

def isolve(a,b,d):
    q, r = divmod(a, b)
    if r == 0:
        return ([0, d/b])
    else:
        sol = isolve(b,r,d)
        u = sol[0]
        v = sol[1]
        return ([v, u - q*v])
    
# Iterative Algorithm (xgcd)
def iterative_egcd(a, b):
    x,y, u,v = 0,1, 1,0
    while a != 0:
        q,r = b//a,b%a; m,n = x-u*q,y-v*q # use x//y for floor "floor division"
        b,a, x,y, u,v = a,r, u,v, m,n
    return b, x, y

def modinv(a, m):
    g, x, y = iterative_egcd(a, m) 
    if g != 1:
        raise ZeroDivisionError, "Inverse does not exist."
    else:
        return x % m



# Chinese Remainder Theorem
# USAGE: icrt(a[], n[])
# Following the form of Sage:
# first parameter: list of a_i, second parameterL list of modulo n_i

def icrt(a, n):
    sum = 0
    prod = reduce(lambda a, b: a*b, n)
 
    for a_i, n_i in zip(a, n):
        p = prod / n_i
        sum += a_i * mul_inv(p, n_i) * p
    return sum % prod

 
def mul_inv(a, b):
    b0 = b
    x0, x1 = 0, 1
    if b == 1: return 1
    while (a > 1):
        if b == 0: return 0
        q = a / b
        a, b = b, a%b
        x0, x1 = x1 - q * x0, x0
    if x1 < 0: x1 += b0
    return x1


# The MotherLoad for Number Theory: see Brian Gladman
# 173.254.28.24/~brgladma/number_theory.py

def gcd(a, *r):
  '''
  Greatest Common Divisor of a sequence of values

  >>> gcd(1701, 13979)
  7

  >>> gcd(117, -17883411)
  39

  >>> gcd(3549, 70161,  336882, 702702)
  273
  '''
  for b in r:
    while b:
      a, b = b, a % b
  return abs(a)


def lcm(a, *r):
  '''
  Least Common Multiple Divisor of a sequence of values

  >>> lcm(1701, 13979)
  3396897

  >>> lcm(117, -17883411)
  53650233

  >>> lcm(3549, 70161,  336882, 702702)
  111426753438
  '''
  for b in r:
    a *= b // gcd(a, b)
  return abs(a)

def xgcd(a, b):
  '''
  Euclid's Extended GCD Algorithm

  >>> xgcd(314159265, 271828186)
  (-18013273, 20818432, 7)
  '''
  u, u1 = 1, 0
  v, v1 = 0, 1
  while b:
    q, r = divmod(a,  b)
    a, b = b, r
    u, u1 = u1, u - q * u1
    v, v1 = v1, v - q * v1
  return (u, v, a) if a > 0 else (-u, -v, -a)

def cong(a, p, n):
  '''
  Solve the congruence a * x == p (mod n)

  >>> cong(13, 17, 19)
  13

  >>> cong(17, 19, 23)
  16

  >>> cong(14, 6, 91) is None
  True
  '''
  g = gcd(a, n)
  if p % g > 0:
    return None
  return (xgcd(a, n)[0] * (p // g)) % (n // g)

def inv(a, n):
  '''
  Find the modular inverse of a (mod n)

  >>> inv(271828, 314159)
  898

  >>> inv(314159, 271828)
  271051
  '''
  return None if gcd(a, n) > 1 else xgcd(a, n)[0] % n


def crt(a, m):
  '''
  The Chinese Remainder Theorem (CRT)

  Solve the equations x = a[i] mod m[i] for x

  >>> crt((2, 3, 5, 7), (97, 101, 103, 107))
  96747802
  '''
  def _crt(a, b, m, n):
    d = gcd(m, n)
    if (a - b) % d:
      return None
    x, y, z = m // d, n // d, (m * n) // d
    p, q, r = xgcd(x, y)
    return (b * p * x + a * q * y) % z

  if len(a) == len(m):
    x, mm = a[0], m[0]
    for i in range(1, len(a)):
      x = _crt(a[i], x, m[i], mm)
      if x is None:
        break
      mm = lcm(m[i], mm)
  else:
    raise IndexError
  return x

# Curiosities
# These I imported functions from numpy to make work
# They work right, but I want to research map, reduce, sum, subtract, etc

from numpy import product, subtract, add
from operator import mod

def eea_map(a,b):
    """Extended Euclidean Algorithm for GCD

    This EEA returns a 3-tuple, the 0th element being the gcd(a,b),
    the next two being s,t such that gcd(a,b) = s*a + t*b, i.e.
    two integers which bring the initial a,b within the gcd of
    one another.

    Example:

    >>> eea_map(17,25)            # gcd(a,b) = 1
    [1, 3, -2]
    >>> 3*17 - 2*25           # s = 3, t = -2
    1

    >>> eea_map(29834,8282)       # gcd(a,b) = 2
    [2, -88, 317]
    >>> -88*29834 + 317*8282  # s = -88, t = 317
    2

    With EEA, you can also find the inverse of x mod n, provided
    gcd(x,n)=1.  The inverse is that number which, multiplying
    x, gives a remainder of 1 when divided by n.

    Example:  find x such that (x * 7) mod 4 = 1.

    Computing:

    >>> inverse(7,4)
    3

    i.e. (3*7) mod 4 = 21 mod 4 = integer remainder of 21/4 = 1.
    """
    v1 = [a,1,0]
    v2 = [b,0,1]
    while v2[0]<>0:
       p = v1[0]//v2[0] # floor division
       v2, v1 = map(subtract,v1,[p*vi for vi in v2]), v2
    return v1

def inverse(m,k):
     """
     Return b such that b*m mod k = 1, or 0 if no solution
     """
     v = eea_map(m,k)
     return (v[0]==1)*(v[1] % k)

def xcrt(a, m):
     """
     Chinese Remainder Theorem:
     m = list of pairwise relatively prime integers
     a = remainders when x is divided by m
     (ai is 'each in a', mi 'each in m')

     The solution for x modulo M (M = product of m) will be:
     x = a1*M1*y1 + a2*M2*y2 + ... + ar*Mr*yr (mod M),
     where Mi = M/mi and yi = (Mi)^-1 (mod mi) for 1 <= i <= r.

     The Chinese Remainder Theorem states that if I give you
     a smattering of divisors, all coprime to each other, and
     tell you what the remainders for each of them is, then
     you can tell me a number which meets my specifications.

     Example:  Your divisors are 7,11 and 15.  The respective
     remainders are 2, 3 and 0.  What's a number that works?

     Computing:

     >>> xcrt([2,3,0], [7,11,15])
     135

     Check:

     >>> 135 % 7
     2
     >>> 135 % 11
     3
     >>> 135 % 15
     0 

     Works.
     """

     M  = product(m)        # multiply m together
     Ms = [M/mi for mi in m]   # list of all M/mi
     ys = [inverse(Mi, mi) for Mi,mi in zip(Ms,m)] # uses inverse,eea
     return sum([ai*Mi*yi for ai,Mi,yi in zip(a,Ms,ys)]) % M

# def prod( iterable ):
#    p = 1
#    for n in iterable:
#        p *= n
#        return p
