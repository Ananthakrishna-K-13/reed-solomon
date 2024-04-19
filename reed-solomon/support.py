import gmpy2
import random

# Converts a string to a mpz bigint object
def bigint(s):
    return gmpy2.mpz(s)

# Generates and returns a random integer in the range [a,b)
def rnd(a, b):
    return bigint(a + random.random() * (b-a))

# Returns the remainder when a is divided by b
def mod(a, b):
    ans = a - b * (a // b)
    if ans < 0 and b > 0:
        ans += b
    if ans < 0 and b < 0:
        ans -= b
    return ans

# Returns a^b
def power(a, b):
    ans = 1
    while b > 0:
        if b & 1:
            ans = ans * a
        a = a * a
        b //= 2
    return ans

# Returns (a^b) mod n
def powermod(a, b, n):
    a = mod(a, n)
    ans = 1
    while b > 0:
        if b & 1:
            ans = mod(ans * a , n)
        a = mod(a * a , n)
        b //= 2
    return ans

# Returns the list of r,s,t in the egcd algorithm. Last element of r is the gcd
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

# Modular inverse of x mod N if gcd(x,N) is 1
def modinv(x, N):
    gcd, ans = egcd(x, N)[0][-1], egcd(x,N)[1][-1]
    if gcd == 1:
        return mod(ans, N)
    raise ValueError("Inverse doesn't exist")

# Chinese remainder theorem
def crt(moduli, remainders):
    N, ans = 1, 0
    for n in moduli:
        N *= n
    for mod, rem in zip(moduli, remainders):
        Ni = N // mod
        Mi = modinv(Ni, mod)
        ans += rem * Ni * Mi
    return ans % N

# Miller - Rabin test for checking primality
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

# Returns n random primes between a and b
def getkprimes(a, b, n):
    primes = []
    while len(primes) < n:
        p = rnd(a, b)
        if miller_rabin(p, 10) and p not in primes:
            primes.append(p)
    return primes

# Returns n random numbers between a and b
def getknumbers(a, b, n):
    ans = []
    for _ in range(n):
        ans.append(bigint(a + random.random() * (b - a)))
    return ans