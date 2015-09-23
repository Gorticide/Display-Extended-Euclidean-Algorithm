# Michael Hentrich is altering this script so as to create pedagical
# version which will display each step of the process.
# euclid(a,b)
# ext_euclid(a,b)
# xed(a,b) for displaying steps of the process
# Date: 9/22/2015
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
