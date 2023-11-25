# TODO
## I suspect the field K.<a> = NumberField(x^4+3x^2-7x + 4) does not satisfy the HNP.
## specifically for 5^2: all local extensions are either 1 or 2 (sometimes 3, but doesnt matter for product formula)
## but im not sure whether there exists an element of norm 5^2.

# aux functions
def flip(A):
    '''
        Flips a matrix A
    '''
    return A[::-1,:].T[::-1,:].T

def col_hermite_form(A):
    '''
        Returns the COLUMN hermite form of given matrix A
        together with its transformation matrix.
        hermite_form in Sage is not sufficient as it uses another
        convention for the hermite form not compatible with Cohen.
    '''
    H, U = flip(A).T.hermite_form(transformation=True)
    return flip(H.T), flip(U.T)

def better_valid_check(a,b):
    '''
        Check whether the field Q(sqrt(a), sqrt(b)) satisfies the Hasse norm principle (HNP).
        Per Tate (Cassels 1967) it is known that such a field satisfies the HNP if and only if
        all its decomposition groups are cyclic.
    '''
    L.<a,b> = NumberField([x^2-a, x^2-b])
    K1.<a1> = L.absolute_field()
    ram_primes = K1.discriminant().prime_divisors()
    # unramified primes are split in at least one subextension.
    # hence their decomposition group is always cyclic
    for p in ram_primes + [2]:
        e, f, _ = K1.decomposition_type(p)[0]
        if e * f == 4:
            return 0
    # also check 2
    return 1

def right_mul(X, A):
    '''
        Compute the multiplicative vector multiplication of Y = XA
        Input:
            X: k-tuple;
            A: (k x l)-matrix with integer values.
        Output:
            Y: l-tuple such that
                Y_j = prod_{i} X_i^{a_{ij}}
    '''
    assert(len(X) == A.nrows())
    return [prod(X[ind]^a for ind, a in enumerate(col)) for col in A.columns()]

# ====================================================

def generate_tuples(N):
    Ps = prime_range(3, N)
    valids = []
    for p, q in Combinations(Ps, k=2):
        if (p % 4 == 1 and q % 4 == 1):
            if legendre_symbol(p, q) == 1:
                valids.append((p,q))
    return valids

def example_biquadratic():
    '''
        Example 7.2.8 of thesis.
    '''
    x = polygen(QQ)
    L.<a> = NumberField(x^4-2*x^3-87*x^2+88*x-22)
    K = L.subfields()[0][0] # Q
    UL = L.unit_group()
    p1, p2 = [I for I in L.ideal(11).prime_factors()]
    S = [p1, p2]


    P = Matrix(ZZ, [I.ideal_class_log() for I in S]).T
    A = block_matrix([P, Matrix(ZZ, [2])], ncols=2)
    _, U = col_hermite_form(A)
    U.subdivide(len(S), len(S))
    U1 = U.subdivision(0,0)
    H, _ = col_hermite_form(U1)
    X = right_mul(S, H)
    t4, t5 = [el.gens_reduced()[0] for el in X]
    t0, t1, t2, t3 = UL.gens_values()
    gens = [t0, t1, t2, t3, t4, t5]
    return gens

def example_relative(a0, a1):
    '''
        Solves N_{Q(zeta_9)/Q(zeta_3)}(x) = a_0 + a_1 zeta_3.
        If it has a solution, returns a random one based on the output
        of the coordinate_vector() method of an integral lattice.
        If no solution, returns 0.
    '''
    # setup
    K.<z3> = CyclotomicField(3)
    x = polygen(K)
    L.<z9> = K.extension(x^3-z3)
    mK = K.ideal(a0 + a1 * z3)
    mL = L.ideal(mK)
    SK = mK.prime_factors()
    SL = mL.prime_factors()
    SUK = K.S_unit_group(S=tuple(SK))
    SUL = L.S_unit_group(S=tuple(SL))

    # calculating matrix of exponents
    ts = SUL.gens_values()
    ss = SUK.gens_values()
    B = Matrix(ZZ, len(ts), len(ss))
    for i, ti in enumerate(ts):
        B[i, :] = vector(SUK.log(ti.relative_norm()))
    A = vector(SUK.log(a0 + a1 * z3))
    # attempt to calculate integer solution
    B = B.insert_row(B.nrows(), [6] + [0] * (B.ncols() - 1))
    H, U = col_hermite_form(B.T)
    k = sum(1 for r in H.columns() if r == 0)
    H = H.delete_columns(range(k))
    U.subdivide(len(ts), k)
    U1 = U.subdivision(0, 0)
    U2 = U.subdivision(0, 1)
    try:
        Z2 = H.solve_right(A)
        sol = U2 * Z2
        for n in sol:
            if n not in ZZ:
                return "No solutions exists"
        return sol, U1, ts
    except ValueError:
        return "No solution exists"
