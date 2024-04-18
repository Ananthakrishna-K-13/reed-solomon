import gmpy2

rs = gmpy2.random_state(hash(gmpy2.random_state()))

# Converts a string to a mpz bigint object
def bigint(s):
    return gmpy2.mpz(s)

# Generates and returns a random integer in the range [a,b)
def rnd(a, b):
    return gmpy2.mpz_random(rs, b - a) + a

# Returns the remainder when a is divided by b
def mod(a, b):
    ans = a - b*(a//b)
    if ans < 0 and b > 0:
        ans += b
    if ans < 0 and b < 0:
        ans -= b
    return ans

def pown(a,b,m):
    a = mod(a,m)
    ans = 1
    while b > 0:
        if b & 1:
            ans = mod(ans * a , m)
        a = mod(a * a , m)
        b //=2
    return ans

# Returns the GCD d of 2 positive integers a and b, and integers s and t such that as + bt = d
def egcd(a, b):
    r, r_dash, e = a, b, 0
    while mod(r, 2) == 0 and mod(r_dash, 2) == 0:
        r, r_dash, e = r//2, r_dash//2, e + 1
    a_dash, b_dash, s, t, s_dash ,t_dash = r, r_dash, 1, 0, 0, 1
    while r_dash != 0:
        while mod(r, 2) == 0:
            r //= 2
            if mod(s, 2) == 0 and mod(t, 2) == 0:
                s, t = s//2, t//2
            else:
                s, t = (s + b_dash)//2, (t - a_dash)//2
        while mod(r_dash, 2) == 0:
            r_dash = r_dash//2
            if mod(s_dash, 2) == 0 and mod(t_dash, 2) == 0:
                s_dash, t_dash = s_dash//2, t_dash//2
            else:
                s_dash, t_dash = (s_dash + b_dash)//2, (t_dash - a_dash)//2
        if r_dash < r:
            r, s, t, r_dash, s_dash, t_dash = r_dash, s_dash, t_dash, r, s, t
        r_dash, s_dash, t_dash = r_dash-r, s_dash-s, t_dash-t
        d = r << e
    return [d,s,t]

# Implementation of CRT using EGCD for the inverse
def crt(n_list, a_list):
    N, x = 1, 0
    N_list, N_bar_list = [], []
    for n in n_list:
        N = N*n
    for n in n_list:
        Ni = N//n
        N_list.append(Ni)
        temp = egcd(Ni, n)[1]
        if(temp < 0) :
            temp = mod(egcd(Ni, n)[1], n)
        N_bar_list.append(temp)
    for i in range(len(n_list)):
        x += a_list[i]*n_list[i]*N_bar_list[i]
    return mod(x, N)

# Miller-Rabin test for checking primality
def miller_rabin(n, rounds):
    if n in [2, 3]:
        return True
    if n==1 or mod(n,2) == 0:
        return False
    t = n - 1
    h = 1
    while mod(t,2) == 0:
        t = t // 2
        h += 1
    for _ in range(rounds):
        a = rnd(2, n-1)
        x = pown(a, t, n)
        for _ in range(h):
            y = pown(x, 2, n)
            if y == 1 and x != 1 and x != n-1:
                return False
            x = y
        if y != 1:
            return False
    return True