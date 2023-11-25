# Attempts to solve the norm form
#   N_{K/Q}(x + y a)
# where a is a generator of K over Q. Due to the cost of casting
# to a Galois closure, it is recommended to only use this program
# for extensions whose Galois closure is of relatively small degree
# i.e. smaller than 7.
prec = 10e-15
R.<x> = QQ[]
f = x^4 - 5
K.<a> = f.root_field()
L.<b> = K.galois_closure()
s,t = K.signature()
n = s + 2 * t
psi = L.complex_embeddings()[0] # fix an embedding of L
galois_embeds = list(K.embeddings(L))
embeds = list(K.real_embeddings())
complex_embeds = K.complex_embeddings()
conjugates = []

# Put the embeddings in a canonical ordering:
#   - first s embeddings are real;
#   - after, t non-conjugated complex embeddings followed by their
#       conjugation
for i, phi in enumerate(complex_embeds):
    if not phi(a) in [s(a) for s in embeds]:
        for phi_conjugate in complex_embeds[i:]:
            if phi_conjugate(a) == phi(a).conjugate():
                embeds.append(phi)
                conjugates.append(phi_conjugate)
embeds = embeds + conjugates
del(complex_embeds, conjugates)

# swap the galois embedding accordingly such that the i-th index
# corresponds to embeds[i] after composition by psi.
for i, phi in enumerate(embeds):
    for j, rho in enumerate(galois_embeds):
        a_approx = rho.post_compose(psi)(a)
        if abs(a_approx - phi(a)) < prec:
            galois_embeds[i], galois_embeds[j] = galois_embeds[j], galois_embeds[i]

OK = K.unit_group()
etas = OK.fundamental_units() # output is random
indices = Combinations(range(s + t), k=s+t-1)
c4 = -Infinity
for index in indices:
    U = Matrix([[log(abs(embeds[i](etas[j]))) for j in range(s+t-1)] for i in index])
    if (U^(-1)).norm(Infinity) > c4:
        c4 = (U^(-1)).norm(Infinity)

c1 = Infinity
c2 = -Infinity
for i, j, k in Permutations(range(n), 3):
    temp1 = abs(embeds[i](a) - embeds[j](a))
    if temp1 < c1:
        c1 = temp1
    temp2 = abs((embeds[j](a) - embeds[k](a)) / (embeds[k](a) - embeds[i](a)))
    if temp2 > c2:
        c2 = temp2
c1 = 2 * c1^(-1)
c3 = c1 * c2
c5 = floor(c4 / (n - 1) * 10^2) / 10^2

# Hardcoded for now, should loop over all values of i
#   and choose j and k most optimal
i, j, k = (0, 2, 1)

# Compute the constant in Baker's theorem
d = L.absolute_degree()
c8 = 18 * factorial(n + 1) * n^(n + 1) * (32*d)^(n+2) * log(2*n*d) * max(1/d, abs(log(-1))/d)
alpha2 = (galois_embeds[i](a) - galois_embeds[j](a)) / (galois_embeds[k](a) - galois_embeds[i](a))
c8 *= max(alpha2.global_height(), 1/d, abs(log(psi(alpha2))) / d)
linear_form = []
for eta in etas:
    temp = galois_embeds[k](eta) / galois_embeds[j](eta)
    linear_form.append(temp)
    c8 *= max(temp.global_height(), 1/d, abs(log(psi(temp))) / d)

X0 = 2 / c5 * (log(2*c3) + c8 * log((s+t)/2) + c8 * log(c8/c5))
# assume the linear form is not real or completely imaginary!
imag_condition = all(imag(log(psi(alph))) == 0 for alph in linear_form)
if imag_condition:
    expon = log_b(X0^(s + t), 10)
    rows = 1
else:
    expon = log_b(X0^((s + t) / 2), 10)
    rows = 2
C = 10^(ceil(expon / 10^floor(log_b(expon, 10))) * 10^floor(log_b(expon, 10)))
new_rows = Matrix(ZZ, rows, s + t)
new_rows[-1, -1] = round(2 * C * pi)
for i, alph in enumerate(linear_form):
    new_rows[0, i] = round(C * real(log(psi(alph))))
    if not imag_condition:
        new_rows[1, i] = round(C * imag(log(psi(alph))))
M = matrix.identity(s + t)
M[s+t-rows:] = new_rows

y = Matrix(ZZ, s + t, 1)
if imag_condition & (imag(psi(alpha2)) == 0):
    y[-1] = -round(C * psi(alpha2))
else:
    y[-2] = -round(C * real(psi(alpha2)))
    y[-1] = -round(C * imag(psi(alpha2)))
ML = M.transpose().LLL().transpose()
MLgram = ML.transpose().gram_schmidt()[0].transpose()
b1 = ML.column(0)
c1 = max(b1.norm()^2 / bi.norm()^2 for bi in MLgram.columns())
sigma = (ML^(-1) * y).column(0)
for i in range(sigma.length() - 1, -1, -1):
    if sigma[i] != 0:
        fracsigma = sigma[i].n().frac()
        if fracsigma < 0:
            # convention of sage: fractional part are negative when the number is
            fracsigma = 1 + fracsigma
        break
c42 = fracsigma * b1.norm()^2
S = (s + t - 2) * X0^2
T = (1 + (s + t) * X0) / sqrt(2)
if c42 > T^2 + S:
    new_bound = (log(C * 2 * c3) - log(sqrt(c42 - S) - T)) / c5
    print(f"New bound: X0 <= {new_bound.n()}")
    print(f"Previous bound: X0 <= {X0.n()}")
else:
    print("c4^2 does not satisfy the desired inequality "
          "try again with a different value of C")
    print(f"c42={c42.n()}")
    print(f"T^2+S={(T^2+S).n()}")
