import gmpy2
import random
import math
from support import *

# M = int(input())
# mu = int(input())
# k = int(input())

M=1000000000000000000
mu = 0.4
k = 12

primes = []
currupted_residues = mu*k



def GlobalSetup(mu, M):
    global primes
    primes = getkprimes(math.sqrt(M), M,k)


def ReedSolomonSend(a):
    tosend = [ a%pi for pi in primes]
    b = Transmit(tosend)
    # print(tosend)
    return b


def Transmit(tosend):
    I = random.sample(range(k),random.randint(0,int(currupted_residues))) 
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

def ReedSolomonReceive(b):
    tempb = crt(primes,b)
    # print(primes)
    # print(b)
    # print(tempb)

    primes.sort(reverse=True)
    P = 1
    n = 1
    for i in range(k):
        if(i<currupted_residues):
            P=P*primes[i]
        n = n * primes[i]
    assert (n > 2*M*P*P)

    rstar = M*P

    lists = extgcd(n,tempb)
    # print(lists[0])
    # print(rstar)
    # print(lists[1])
    # print(lists[2])
    # print(rstar)

    ind = 0
    while(lists[0][ind]>rstar):
        ind = ind+1

    rdash,sdash,tdash = lists[0][ind],lists[1][ind],lists[2][ind]
    # mod(tdash,n)

    # n,tempb,rstr,P  ==> rdash,sdash,t_dash


    if(rdash%tdash == 0):
        return(rdash//tdash)
    else:
        return("Error")
    






GlobalSetup(mu,M)
print(ReedSolomonReceive(ReedSolomonSend(34664363434737)))






