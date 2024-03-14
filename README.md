# FGLM from scratch

Implementing the FGLM (Faugère, Gianni, Lazard, and Mora) algorithm for converting Gröbner bases of zero-dimensional ideals with respect to one term/monomial order to another.

Work in progress. To run unit tests:

    ./test.py

Example usage:

    $ python
    >>> from polynomials import *
    >>> str(Polynomial.parse("3x^2 - x + 1", lex))
    '3x^2 - x + 1'
    >>> str(Polynomial.parse("3x^2 - x + 1", lex))
    '3x^2 - x + 1'
    >>> str(Polynomial.parse("x^2 + 2x + 1", lex) / Polynomial.parse("x + 1", lex))
    'x + 1'
    >>> from groebner import *
    >>> for p in quotient_monomials([Polynomial.parse(p, lex) for p in "x^2 + y^2 - 1, x^2 + z^2 - 1".split(",")], max_degree=2):
    ...     print(p)
    ... 
    x
    z
    y
    x^2 + y^2 - 1
    -y^2 + 1
    xz
    xy
    z^2
    yz
    y^2
    >>> 
