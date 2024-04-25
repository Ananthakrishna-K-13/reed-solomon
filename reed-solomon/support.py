import gmpy2
import random

'''
Support functions for reed_solomon.py

bigint(s) -> Converts a string to a mpz bigint object
rnd(a,b) -> Returns a random integer in range [a,b)
mod(a,b) -> Returns a % b, always in range 0 to b-1
power(a,b) -> Returns a**b
powermod(a,b,n) -> Returns a**b % n
egcd(a,b) -> Returns a list of 3 lists r,s,t, such that as + bt = r for each index. The last entry of r is the gcd
modinv(x,N) -> Returns inverse of x modulo N if their gcd is 1, using egcd
crt(moduli,remainders) -> Returns solution to linear congruences x = moduli mod remainder
miller_rabin(n,rounds) -> Returns if a number is prime or not by doing n rounds of Miller-Rabin test
'''

def bigint(s):
    return gmpy2.mpz(s)

def rnd(a, b):
    return bigint(a + random.random() * (b-a))

def mod(a, b):
    ans = a - b * (a // b)
    if ans < 0 and b > 0:
        ans += b
    if ans < 0 and b < 0:
        ans -= b
    return ans

def power(a, b):
    ans = 1
    while b > 0:
        if b & 1:
            ans = ans * a
        a = a * a
        b //= 2
    return ans

def powermod(a, b, n):
    a = mod(a, n)
    ans = 1
    while b > 0:
        if b & 1:
            ans = mod(ans * a , n)
        a = mod(a * a , n)
        b //= 2
    return ans

def egcd(a,b):
    r, r_dash, e = a, b, 0
    r_list, s_list, t_list = [], [], []

    while mod(r, 2) == 0 and mod(r_dash, 2) == 0:
        r, r_dash, e = r // 2, r_dash // 2, e + 1
    a_dash, b_dash, s, t, s_dash ,t_dash = r, r_dash, 1, 0, 0, 1

    while r_dash != 0:
        while mod(r, 2) == 0:
            r //= 2
            if mod(s, 2) == 0 and mod(t, 2) == 0:
                s, t = s // 2, t // 2
            else:
                s, t = (s + b_dash) // 2, (t - a_dash) // 2
            r_list.append(r * power(2, e))
            s_list.append(s)
            t_list.append(t)

        while mod(r_dash, 2) == 0:
            r_dash //= 2
            if mod(s_dash, 2) == 0 and mod(t_dash, 2) == 0:
                s_dash, t_dash = s_dash // 2, t_dash // 2
            else:
                s_dash, t_dash = (s_dash + b_dash) // 2, (t_dash - a_dash) // 2
            r_list.append(r_dash * power(2, e))
            s_list.append(s_dash)
            t_list.append(t_dash)

        if r_dash < r:
            r, s, t, r_dash, s_dash, t_dash = r_dash, s_dash, t_dash, r, s, t
        r_dash, s_dash, t_dash = r_dash - r, s_dash - s, t_dash - t

    return [r_list, s_list, t_list]

def modinv(x, N):
    gcd, ans = egcd(x, N)[0][-1], egcd(x,N)[1][-1]
    if gcd == 1:
        return mod(ans, N)
    raise ValueError("Inverse doesn't exist")

# def crt(moduli, remainders):
#     N, ans = 1, 0
#     for n in moduli:
#         N *= n
#     for mod, rem in zip(moduli, remainders):
#         Ni = N // mod
#         Mi = modinv(Ni, mod)
#         ans += rem * Ni * Mi
#     return ans % N

def crt(primes_list,a_list):
    n=gmpy2.mpz(1)
    for i in primes_list:
        n=n*i
    k=len(primes_list)
    nstar=[gmpy2.mpz(1) for i in range(k)]
    e=[gmpy2.mpz(1) for i in range(k)]
    for i in range(k):
        nstar[i]=n//primes_list[i]
        b=mod(nstar[i],primes_list[i])
        t=modinv(b,primes_list[i])
        e[i]=nstar[i]*t
    ans=gmpy2.mpz(0)
    for i in range(k):
        ans=mod(ans+a_list[i]*e[i],n)
    return ans

def miller_rabin(n, rounds):
    if n in [2, 3]:
        return True
    if n == 1 or mod(n, 2) == 0:
        return False
    t, h = n - 1, 1
    while mod(t, 2) == 0:
        t = t // 2
        h += 1
    for _ in range(rounds):
        a = rnd(2, n-1)
        x = powermod(a, t, n)
        for _ in range(h):
            y = powermod(x, 2, n)
            if y == 1 and x != 1 and x != n-1:
                return False
            x = y
        if y != 1:
            return False
    return True

def getkprimes(n):
    primes = []
    s = 7993632487
    while len(primes) < n:
        if miller_rabin(s,10):
            primes.append(s)
        s += 2
    return primes

def getknumbers(a, b, n):
    ans = []
    for _ in range(n):
        ans.append(bigint(a + random.random() * (b - a)))
    return ans

def extgcd(a,b):
    r,rdash,s,sdash,t,tdash = a,b,1,0,0,1
    rlist,slist,tlist = [r],[s],[t]
    while(rdash!=0):
        q,rddash = r//rdash,mod(r,rdash)
        r,s,t,rdash,sdash,tdash = rdash,sdash,tdash,rddash,s-sdash*q,t-tdash*q
        rlist.append(r)
        slist.append(s)
        tlist.append(t)

    return [rlist,slist,tlist]