import gmpy2
import random
import math
from support import *


M=gmpy2.mpz('1'+'0'*1000)
mu = random.random()
# mu = 0.43
k = 1000

primes = []
currupted_residues = mu*k


def GlobalSetup(mu, M):
    global primes
    primes = getkprimes(k)


def ReedSolomonSend(a):
    tosend = [ mod(a,pi) for pi in primes]
    b = Transmit(tosend)
    return b


def Transmit(tosend):
    I = random.sample(range(k),random.randint(0,int(currupted_residues))) 
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

    primes.sort(reverse=True)
    P = 1
    n = 1
    for i in range(k):
        if(i<currupted_residues):
            P=P*primes[i]
        n = n * primes[i]

    if(n <= 2*M*P*P):
        return(-1)

    rstar = M*P
    lists = extgcd(n,tempb)

    ind = 0
    while(lists[0][ind]>rstar):
        ind = ind+1

    rdash,sdash,tdash = lists[0][ind],lists[1][ind],lists[2][ind]

    if(rdash%tdash == 0):
        return(rdash//tdash)
    else:
        return(-1)
    
def main():
    GlobalSetup(mu,M)
    a = random.randint(0,M)
    print("Message transmitted:")
    print(a)
    stat = ReedSolomonReceive(ReedSolomonSend(a))
    # print(stat)
    if(stat == -1):
        print(mu)
        print("Message could not be reconstructed")
    else:
        print("Message recieved:")
        print(stat)


main()













