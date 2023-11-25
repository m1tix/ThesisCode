from sage.modules.free_module_integer import IntegerLattice

def squarefree_gen(n):
    '''
        Return all squarefree positive integers <= n by sieving
    '''
    nums = [1] * (n + 1)
    # we consider 1 not squarefree :)
    nums[0] = 0
    nums[1] = 0
    for i in range(2, ceil(sqrt(n))):
        if nums[i]:
            for m in range(i**2, n + 1, i ** 2):
                nums[m] = 0
    return [i for i, b in enumerate(nums) if b]


def soluble_locally_everywhereD(D, q_max=10000):
    '''
        Given D, attempts to find a prime q such that the equation
            x^2 - Dy^2 = q
        is soluble modulo all primes p. User can increase the search range of q
        by specifying q_max.
    '''
    if D % 4 == 1:
        return None
    D_primes = prime_factors(D)
    if D_primes[0] == 2:
        del D_primes[0]
    if D % 8 == 2:
        congr_cond = [1, 7]
    elif D % 8 == 6:
        congr_cond = [1, 3]
    else:
        congr_cond = [1, 5]
    for q in prime_range(3, q_max):
        if not q % 8 in congr_cond:
            continue
        for p in D_primes:
            if legendre_symbol(q, p) != 1:
                break
        else:
            return q
    return None

def generate_random(N):
    '''
        Generate a 'random' quadratic field K = Q(sqrt(D)) with D <= N, D != 1 mod 4
        and a prime number q such that
            i) N_{K/Q}(x) = q is non-soluble in O_K;
            ii) N_{K/Q}(x) = q is soluble modulo all rational primes.
        Since HNT applies, also returns an element of O_K of norm q.
    '''
    # Find D <= N and q that satisfy above conditions
    sqfrs = squarefree_gen(N)
    sqfrs = [s for s in sqfrs if s % 4 != 1]
    valid = False
    while not valid:
        D = sample(sqfrs, k=1)[0]
        q = soluble_locally_everywhereD(D)
        L.<a> = QuadraticField(D)
        temp_I = L.ideal(q).prime_factors()[0]
        if not temp_I.is_principal():
            valid = True
    # Compute necessary invariants of L
    K = L.subfields()[0][0]
    clL = L.class_group()
    S = [q]
    T = L.ideal(q).prime_factors()
    for I in clL.gens_ideals():
        temp = prime_factors(I.norm())
        S += temp
        for p in temp:
            T += L.ideal(p).prime_factors()
    ULT = L.S_unit_group(S=tuple(T))
    ts = ULT.gens_values()
    UKT = K.S_unit_group(S=tuple(S))
    ss = UKT.gens_values()
    
    # Use Algorithm 7.2.6 to find a solution (existence is guarenteed)
    B = matrix(ZZ, len(ts), len(ss))
    for i, ti in enumerate(ts):
        B[i, :] = vector(UKT.log(ti.norm()))
    A = vector(UKT.log(q))
    B = B.insert_row(B.nrows(), [2] + [0] * (B.ncols() - 1))
    BLLL, U = B.LLL(transformation=True)
    nz = sum(1 for r in BLLL.rows() if r == 0)
    BLLL = BLLL.delete_rows(range(nz))
    U = U.delete_rows(range(nz))
    lat_B = IntegerLattice(BLLL, lll_reduce=False)
    v = lat_B.coordinate_vector(A) * U
    x = prod(ti^v[i] for i, ti in enumerate(ts))
    return (D, q, x)
