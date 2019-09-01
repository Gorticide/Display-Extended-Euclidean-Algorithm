# Display-Extended-Euclidean-Algorithm
This prints out each step of the extended Euclidean algorithm
 Michael Hentrich is altering this script so as to create pedagogical
 version which will display each step of the process.
 Date: September 2015

 xed(a,b) for displaying steps of the process
 xlgcd(a,b) - like xed but includes detailed summary
 isolve(a,b,d) - solves Diophntine equations using euclid
 iterative_egcd(a, b)
 modinv(a, m)
 icrt(a, n) - Chinese Remainder Theorem
 mul_inv(a, b) - multiplicative modular inverse

 The MotherLoad for Number Theory: see Brian Gladman
 173.254.28.24/~brgladma/number_theory.py
 crt, gcd, xgcd, lcm, cong, inv
 
 Kirby Urner for a cooler version of the Extended Euclidean Algorithm, with
 related inverse and Chinese Remainder Theorem functions.
 euclid(a,b)
 ext_euclid(a,b)

 Functions at the the end are curosities where I am testing out
 using product, add, subtract from numpy to use with reduce()
 eea_map, xcrt, inverse

 Author: Sam Erickson for test(a,b) and getCoeffic(a,b)
 Date: 7/14/2013

 Program Description: This program gives the integer coefficients x,y to the
 equation ax+by=gcd(a,b) given by the extended Euclidean Algorithm. Note that
 this program does not use the typical extended Euclidean Algorithm to find
 ax+by=gcd(a,b).



General algorithm for program:

1. Get integers a and b from user.

2. Call test(a,b) and assign to get one integer coefficient for the equation
   (1/gcd(a,b))*[a*x+b*y=gcd(a,b)].

3. Call getCoffic(min2,aNew,bNew,GCD) to return a string for of the equation
   a*x+b*y=gcd(a,b) for printing.

