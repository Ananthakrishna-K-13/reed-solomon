from support import *
import random
from sympy.ntheory.modular import solve_congruence
from sympy import isprime

def test_egcd():
    for _ in range(10):
        a = int(10**10 + random.random()*(10**11))
        b = int(10**10 + random.random()*(10**11))
        ans = egcd(a,b)
        r,s,t = ans[0],ans[1],ans[2]
        for i in range (len(r)):
            if r[i] != a*s[i] + b*t[i]:
                print("EGCD FAIL")
                return None
        if r[-1]!=gmpy2.gcd(a,b):
            print("EGCD FAIL")
            return None
    print("EGCD PASS")
    return None

def test_crt():
    for _ in range (10):
        a = rnd(10,50)
        primes = getkprimes(1000,100000,a)
        moduli = getknumbers(10,1000,a)
        c1 = crt(primes,moduli)
        c2 = solve_congruence(*zip(moduli,primes))[0]
        if c1 != c2:
            print("CRT FAIL")
            return None
    print("CRT PASS")
    return None

def test_miller():
    for i in range(10):
        ans = getknumbers(100000000,10000000000,100)
        for i in ans:
            if miller_rabin(i,10) != isprime(i):
                print("MILLER FAIL")
                return None
    print("MILLER PASS")
    return None

test_egcd()
test_crt()
test_miller()