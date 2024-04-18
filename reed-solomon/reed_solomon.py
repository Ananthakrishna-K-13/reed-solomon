import gmpy2
import random
import math
from support import *

# M = int(input())
# mu = int(input())
# k = int(input())

M=50
mu = 0.42
k = 6

primes = []
currupted_residues = mu*k



def GlobalSetup(mu, M):
    global primes
    primes = getkprimes(math.sqrt(M), 1000000007,k)


def ReedSolomonSend(a):
    tosend = [ a%pi for pi in primes]
    b = Transmit(tosend)
    print(tosend)
    return b


def Transmit(tosend):
    I = random.sample(range(k),random.randint(0,int(currupted_residues))) # why mu*k ??
    print(I)
    b = []
    for i in range(k):
        if(i in I):
            toappend = tosend[i]
            while(toappend==tosend[i]):
                toappend = random.randint(0,primes[i]-1)
            b.append(toappend)
        else:
            b.append(tosend[i])
    return b

# def ReedSolomonReceive():




GlobalSetup(mu,M)
print(ReedSolomonSend(100000000723566))






