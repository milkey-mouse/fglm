#!/usr/bin/env python3
from collections import Counter
from itertools import combinations_with_replacement, count

from polynomials import Term, Polynomial, Lex, DegLex, DegRevLex


def gcd(a, b):
    """Compute the greatest common divisor of two integers."""
    while b != 0:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Compute the least common multiple of two monomials."""
    if isinstance(a, int) and isinstance(b, int):
        return abs(a * b) // gcd(a, b)
    elif isinstance(a, Term) and isinstance(b, Term):
        if a.ordering != b.ordering:
            raise ValueError(
                f"Orderings differ: {self.ordering} != {other.ordering}"
            )

        def vars(a, b):
            i = j = 0
            while i < len(a.vars) and j < len(b.vars):
                term_a = (var_a, exp_a) = a.vars[i]
                term_b = (var_b, exp_b) = b.vars[j]
                if var_a == var_b:
                    yield (var_a, max(exp_a, exp_b))
                    i += 1
                    j += 1
                elif term_a < term_b:
                    yield a.vars[i]
                    i += 1
                else:
                    yield b.vars[j]
                    j += 1
            yield from a.vars[i:]
            yield from b.vars[j:]

        coefficient = lcm(a.coefficient, b.coefficient)
        return Term(coefficient, tuple(vars(a, b)), a.ordering)
    else:
        raise TypeError(f"Unsupported types: {type(a)}, {type(b)}")


def S(p, q):
    """Compute the S-polynomial for p and q."""
    L = lcm(p.leading_term(), q.leading_term())
    return (L / p.leading_term()) * p - (L / q.leading_term()) * q


# I wish someone had told me Buchberger's algorithm is just Knuth-Bendix...
def buchberger(F):
    """Compute the reduced Gröbner basis of a set of polynomials."""
    G = list(F)

    # find a Gröbner basis by repeatedly adding reduced S-polynomials to G
    pairs = [(p, q) for p in G for q in G if p != q]
    while pairs:
        p, q = pairs.pop()
        _, r = S(p, q).reduce(G)

        if r.terms:
            pairs.extend((g, r) for g in G)
            G.append(r)

    # make G a *reduced* Gröbner basis by reducing all elements wrt each other
    for i, g in enumerate(G):
        _, r = g.reduce(G[:i] + G[i + 1 :])
        G[i] = r

    return G


def quotient_monomials(G, max_degree=None):
    """
    All monomials in the quotient ring K[x_1, ..., x_n] / <G> up to the given
    max_degree. G must be a Gröbner basis. If max_degree is not provided, <G>
    must be a zero-dimensional ideal, and all monomials in the quotient ring
    (i.e. all monomials in K[x_1, ..., x_n] up to reduction by G) are yielded.
    """
    ordering = next(iter(G)).ordering
    assert all(g.ordering == ordering for g in G)

    all_vars = set(
        (var, 1) for g in G for term in g.terms for var, _ in term.vars
    )
    g_degree = max(term.degree() for g in G for term in g.terms)

    seen = set()
    seen_degree = None
    for degree in range(1, max_degree + 1) if max_degree else count(1):
        for vars in combinations_with_replacement(all_vars, degree):
            _, normal_form = Term(1, vars, ordering).reduce(G)
            if normal_form not in seen:
                seen.add(normal_form)
                print(len(seen))
                seen_degree = degree
                yield normal_form
        if not max_degree and seen_degree and degree - seen_degree > g_degree:
            # because any monomial of degree d or higher must be divisible by
            # at least one of the monomials of degree d-n, d-n+1, ..., d-1 with
            # n the number of variables, if no new normal forms were found in
            # the last n degrees, we won't find more & must have found them all
            break


def fglm(G, target_ordering):
    """
    Convert a Gröbner basis G with respect to one monomial ordering to a Gröbner
    basis with respect to another monomial ordering using the FGLM algorithm.
    """
    # Initialize an empty set to store the new Gröbner basis
    G_target = []

    # Initialize the set of standard monomials with the constant monomial 1
    N_target = [Term(1, (), target_ordering)]

    # Iterate over all monomials in the target ordering
    # TODO: I removed TermOrdering.monomials because I think it was wrong,
    # the alternative should maybe involve quotient_monomials when it works
    for monomial in target_ordering.monomials(
        G[0].ordering.vars, max(g.degree() for g in G)
    ):
        # Compute the normal form of the monomial with respect to the input basis G
        _, p = monomial.reduce(G)

        # Extract the coefficients of the normal form with respect to the current set of standard monomials
        coefficients = [
            term.coefficient for term in p.terms if term in N_target
        ]

        # If all coefficients are zero, add the monomial to the set of standard monomials
        if not any(coefficients):
            N_target.append(monomial)
        else:
            # Construct a new polynomial with leading term equal to the current monomial
            g = monomial - p

            # Add the new polynomial to the Gröbner basis
            G_target.append(g)

            # Elimination test: if the leading monomial is a power of the first variable, terminate
            if (
                g.leading_term().vars
                and g.leading_term().vars[0][0] == target_ordering.vars[0]
            ):
                return G_target

    return G_target
