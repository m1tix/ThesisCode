import time


def find_gamma_lambda(s, t):
    '''
        Find appriopriate values of gamma and lambda as defined in
        Pohst p. 342. Note that this function actually returns gamma/m, lam,
        so the return value still needs to be multiplied by m for the actual value of gamma.
    '''
    cache = {
        (0, 0): (0.0, 1.0000000286422652),
        (0, 1): (0.0, 1.0000000286422652),
        (0, 2): (0.672285914742532, 4.283239982605843),
        (0, 3): (1.8238586120577174, 5.451931914985356),
        (0, 4): (3.7803659183937652, 6.095639968780807),
        (0, 5): (7.100618497767398, 6.503256952373093),
        (1, 0): (0.0, 1.0000000286422652),
        (1, 1): (0.6804340401988702, 5.451931835420812),
        (1, 2): (1.8461586087236355, 6.503257075719373),
        (1, 3): (3.826141201818861, 6.990160381827038),
        (1, 4): (7.186720777674836, 7.270851133826022),
        (2, 0): (0.6971589089131371, 8.314635418743713),
        (2, 1): (1.8803482943780323, 8.314635211706706),
        (2, 2): (3.8884088131579393, 8.314635334254142),
        (2, 3): (7.296406522369098, 8.31463533134185),
        (2, 4): (13.080320442645707, 8.314635399161695),
        (3, 0): (1.9392841410811692, 12.117865214745104),
        (3, 1): (3.977962076335489, 10.463180986848297),
        (3, 2): (7.440825137531901, 9.810932268095886),
        (3, 3): (13.318308677294649, 9.462310971537855),
        (4, 0): (4.117501178256238, 14.479077066316238),
        (4, 1): (7.639390967351273, 12.117864956982956),
        (4, 2): (13.623913941787855, 11.06116167885242),
        (4, 3): (23.780106442750412, 10.463180988387565),
        (5, 0): (7.9290446941275725, 16.072848451769744),
        (5, 1): (14.030319172105056, 13.424364592836742),
        (5, 2): (24.393623874611094, 12.117864847603146),
        (6, 0): (14.596245542298588, 17.217192542933255),
        (6, 1): (25.188818310812717, 14.479077066650865),
        (6, 2): (43.164927190565095, 13.020814879054535),
        (7, 0): (26.25859157373109, 18.0774344152063),
        (7, 1): (44.67801290992938, 15.346949366932632),
        (8, 0): (46.660080716738655, 18.74717745410395),
        (8, 1): (78.72783571224137, 16.07284821120593),
        (9, 0): (82.35243085007099, 19.283159040369682),
        (10, 0): (144.80096262226309, 19.721698143341136)
    }
    if s + 2 * t <= 10:
        return cache[(s, t)]
    var('x')
    h(x) = x / (x - 1) - 1 / log(x)
    g(x) = (1 - h(x)) * x ^ h(x) + h(x) * x ^ (h(x) - 1)
    f(x) = (log(x)) ^ (1 - s - t) * g(x) ^ ((s + 2 * t) / 2)
    y_star, l_star = find_local_minimum(f, 1, 5 * (s + 2 * t))
    clear_vars()
    return y_star * (log(l_star)) ^ (s + t - 1) - 1, l_star


def embeddings_of_field(K):
    '''
    Put the embeddings in a canonical ordering:
      - first s embeddings are real;
      - after, t non-conjugated complex embeddings followed by their
          conjugation
    '''
    complex_embeds = K.complex_embeddings()
    embeds = list(K.real_embeddings())
    conjugates = []
    for i, phi in enumerate(complex_embeds):
        if not phi(a) in [s(a) for s in embeds]:
            for phi_conjugate in complex_embeds[i:]:
                if phi_conjugate(a) == phi(a).conjugate():
                    embeds.append(phi)
                    conjugates.append(phi_conjugate)
    return embeds + conjugates


def compute_conj_bounds(K, m):
    '''
        Compute bounds on the conjugates x of norm equation
            |N(x)| = m.
        Ordering is based on embeddings_of_field
    '''
    n = K.degree()
    epsilons = K.unit_group().fundamental_units()
    embeds = embeddings_of_field(K)
    conj_bounds = []
    for i in range(n):
        sum_val = sum(abs(log(abs(embeds[i](ep)))) for ep in epsilons)
        upper_bound = (exp(1 / n * log(m) + 1 / 2 * sum_val)).n()
        lower_bound = (exp(1 / n * log(m) - 1 / 2 * sum_val)).n()
        conj_bounds.append((lower_bound, upper_bound))
    return conj_bounds


def compute_integer_bounds(K, m):
    '''
        Compute bounds on the integers r_j as given in Fincke-Pohst.
    '''
    s, t = K.signature()
    n = s + 2 * t
    gam, lam = find_gamma_lambda(s, t)
    conj_bounds = compute_conj_bounds(K, m)
    int_bounds = []
    for rj, sj in conj_bounds:
        lj = floor(-2 / log(lam) * (log(sj) - log(m) / n))
        uj = ceil(-2 / log(lam) * (log(rj) - log(m) / n))
        int_bounds.append((lj, uj))
    return int_bounds


def compute_valid_forms(K, m):
    '''
        Computes all pairs (r_1, ..., r_n) where n is degree of K such that
            i) r_1 + ... + r_n = 0;
            ii) L_j <= r_j <= U_j, where (L_j, U_j) is the j-th component of 
                the output of compute_integer_bounds;
            iii) if (s,t) is the signature of K then for all 1 <= j <= t
                    r_{s+t+j} - r_{s+j} = 0, 1.
                 Moreover, it is 1 for at most two j's.
    '''
    bounds = compute_integer_bounds(K, m)
    valid_tuples = []
    s, t = K.signature()
    for r in cartesian_product([range(x[0], x[1] + 1) for x in bounds]):
        if sum(r) != 0:
            continue
        count = 0
        for j in range(t):
            val = r[s + t + j] - r[s + j]
            if val == 1:
                count += 1
            elif val != 0:
                break
        else:
            if count < 2:
                valid_tuples.append(r)
    return valid_tuples


def cholesky(A):
    '''
        Returns the UPPER triangular cholesky decomposition of A.
        Assumes that A is positive definite
    '''
    Q = copy(A)
    n = A.nrows()
    for i in range(n):
        for j in range(i + 1, n):
            Q[j, i] = Q[i, j]
            Q[i, j] = Q[i, j] / Q[i, i]
        for k, l in cartesian_product([range(i + 1, n)] * 2):
            Q[k, l] = Q[k, l] - Q[k, i] * Q[i, l]
    R = diagonal_matrix((sqrt(k) for k in Q.diagonal()))
    for i in range(n):
        for j in range(i + 1, n):
            R[i, j] = R[i, i] * Q[i, j]
    return R


def quadratic_supplement(A):
    Q = copy(A)
    n = A.nrows()
    for i in range(n):
        for j in range(i + 1, n):
            Q[j, i] = Q[i, j]
            Q[i, j] = Q[i, j] / Q[i, i]
        for mu in range(i + 1, n):
            for v in range(mu, n):
                Q[mu, v] = Q[mu, v] - Q[mu, i] * Q[i, v]
    for i in range(n):
        for j in range(i):
            Q[i, j] = 0
    return Q


def short_vectors(A, c):
    '''
        Returns all integral vectors v such that Q(v) <= c,
        where Q is the quadratic form defined by A. Only outputs one of
        v, -v and moreover assumes that A is positive definite. Uses Fincke-Pohst
        but should not be used: short_vectors_gp is much faster due to LLL reduction.
    '''
    Q = quadratic_supplement(A)
    n = Q.nrows()
    res = [[]]
    for i in range(n - 1, -1, -1):
        curr = []
        for prev in res:
            x = flatten([[0] * (i + 1), prev])
            si = sum(Q[k, k] * (x[k] + sum(Q[k, j] * x[j]
                                           for j in range(k + 1, n)))**2
                     for k in range(i + 1, n))
            ti = sum(Q[i, j] * x[j] for j in range(i + 1, n))
            ui = floor(sqrt((c - si) / Q[i, i]) - ti)
            li = ceil(-sqrt((c - si) / Q[i, i]) - ti)
            for xi in range(li, ui + 1):
                new_x = [xi] + prev
                mm = -1
                for m in range(len(new_x)):
                    if new_x[m] != 0:
                        mm = m
                if mm == -1:
                    curr.append(new_x)
                else:
                    if new_x[mm] > 0:
                        curr.append(new_x)
        res = curr
    return res


def short_vectors_gp(A, c):
    '''
        Returns all integral vectors v such that Q(v) <= c,
        where Q is the quadratic form defined by A. Only outputs one of
        v, -v and moreover assumes that A is positive definite.
        Uses PARI on the backend.
    '''
    _, _, res = gp.qfminim(A, c, -1, 2)
    return res.sage().columns()


def integral_elements_of_norm(K, m):
    '''
        Returns a complete set of integral elements x mod units of positive norm
        such that 
            |N(x)| = m
        by using class group computations. Should generally not be used,
        K.elements_of_norm(m) is generally much faster.
    '''
    clK = K.class_group()
    factor_m = factor(m)
    ps = [f[0] for f in factor_m]
    es = [f[1] for f in factor_m]
    qs = []
    poss_sols = []
    # determine all solutions in terms of prime ideals
    for i, p in enumerate(ps):
        qp = K.primes_above(p)
        qs.append(qp)
        fs = [q.residue_class_degree() for q in qp]
        # just a test to see if this is faster
        mlp = MixedIntegerLinearProgram()
        w = mlp.new_variable(integer=True, nonnegative=True)
        mlp.add_constraint(sum(fs[j] * w[j] for j in range(len(qp))) == es[i])
        local_sols = mlp.polyhedron().integral_points()
        if local_sols:
            poss_sols.append(local_sols)
        else:
            return []
    # There are solutions in ideals, but these need not be principle.
    # each solution hence must satisfy a certain linear congruence in class group
    dks = clK.gens_orders()  # not elementary_divisors as ordering is different
    aks = []
    for i, p in enumerate(ps):
        curr = []
        for q in qs[i]:
            curr.append(q.ideal_class_log())
        aks.append(curr)
    real_sols = []
    for x in cartesian_product(poss_sols):
        for k, dk in enumerate(dks):
            tot = 0
            for i, row in enumerate(x):
                tot += sum(aks[i][j][k] * xij for j, xij in enumerate(row))
            if tot % dk != 0:
                break
        else:
            # ideal is principal
            ideal_x = prod(q ^ x[i][j] for i, row in enumerate(qs)
                           for j, q in enumerate(row))
            sol = ideal_x.gens_reduced()[0]
            real_sols.append(sol)
            # TODO: check for sign
    return real_sols


def fincke_pohst(K, m, ineq=False, verbose=True):
    '''
        Compute a complete set of representatives of non-associate integral
        elements x of K such that
            |N(x)| = m
        by means of Fincke-Pohst. If ineq is True, then return all integral elements of absolute norm at most m.
        Note that the outcome is not minimal in the sense that each solution may contain
        a non-trivial combination of a norm one element.
    '''
    # Basic notations
    if verbose:
        start_time = time.time()
    s, t = K.signature()
    n = s + 2 * t
    embeds = embeddings_of_field(K)
    gam, lam = find_gamma_lambda(s, t)
    UK = K.unit_group()
    epsilons = UK.fundamental_units()
    ais = K.integral_basis()
    if ineq:
        norms = dict({k: [] for k in range(1, m + 1)})

    # Compute bounds on conjugates of x and bounds on the rj's defining quadratic forms
    bounds = []
    for j in range(n):
        sum_val = sum(abs(log(abs(embeds[j](ep)))) for ep in epsilons)
        print(sum_val)
        lj = floor(-sum_val / log(lam))
        uj = ceil(sum_val / log(lam))
        bounds.append((lj, uj))
    rs = []
    # Compute all tuples subject to the bound above (same as integer_bounds)
    for r in cartesian_product([range(x[0], x[1] + 1) for x in bounds]):
        if sum(r) != 0:
            continue
        count = 0
        for j in range(t):
            val = r[s + t + j] - r[s + j]
            if val == 1:
                count += 1
            elif val != 0:
                break
        else:
            if count < 2:
                rs.append(r)
    quad_bound = ceil(n * m ^ (2 / n) * (1 + gam) ^ (2 / n))
    if verbose:
        inter_time = time.time()
        print(f"Initializiation took: {inter_time - start_time}")
        print(f"Need to check {len(rs)} quadratic forms...")
        print(f"Bound for quadratic form: {quad_bound}")
    res = []
    mats = []
    # iterate over all tuples, each defining a different quadratic form
    for k, r in enumerate(rs):
        if verbose and k != 0:
            print(
                f"Currently at matrix {k}, previous took {time.time() - inter_time }."
            )
            inter_time = time.time()
        B = Matrix(RR, n, n)
        for i, j in cartesian_product([range(n)] * 2):
            if j >= s + t:
                B[j, i] = sqrt(lam ^ r[j] + lam ^ r[j - t]) * imag(
                    embeds[j - t](ais[i]))
            elif j >= s:
                B[j, i] = sqrt(lam ^ r[j] + lam ^ r[j + t]) * real(embeds[j](
                    ais[i]))
            else:
                B[j, i] = lam ^ (r[j] / 2) * embeds[j](ais[i])
        # A is the real symmetric positive definite matrix defining quadratic form
        A = B.T * B
        mats.append(A)
        # Using fincke-pohst, compute all vectors v such that v A v <= quad_bound
        # Then iterate over all of them to check for desired norm
        vecs = short_vectors_gp(A, quad_bound)
        for v in vecs:
            x = sum(ais[i] * v[i] for i in range(n))
            curr_norm = abs(x.norm())
            if curr_norm == 0:
                continue
            if (ineq and curr_norm <= m):
                for y in norms[curr_norm]:
                    if x / y in UK:
                        break
                else:
                    norms[curr_norm].append(x)
            elif curr_norm == m:
                for y in res:
                    if x / y in UK:
                        break
                else:
                    res.append(x)
    if verbose:
        print(f"Final matrix took {time.time() - inter_time}")
        print(f"Total time elapsed: {time.time() - start_time}")
    if ineq:
        return norms
    return res
