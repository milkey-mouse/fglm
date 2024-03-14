from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
from functools import total_ordering
from itertools import combinations_with_replacement
from typing import Callable, Optional
import re, string


@dataclass(frozen=True)
class TermOrdering(ABC):
    vars: tuple[str, ...] | str
    nickname: Optional[str] = field(default=None, repr=False)

    @abstractmethod
    def __call__(self, a: "Term", b: "Term") -> int:
        pass


@dataclass(frozen=True)
class Lex(TermOrdering):
    """Lexicographic term ordering."""

    def __call__(self, a: "Term", b: "Term") -> int:
        for (var_a, exp_a), (var_b, exp_b) in zip(a.vars, b.vars):
            idx_a = self.vars.index(var_a)
            idx_b = self.vars.index(var_b)
            if idx_a < idx_b:
                return -1
            elif idx_a > idx_b:
                return 1
            elif exp_a < exp_b:
                return -1
            elif exp_a > exp_b:
                return 1
        if len(a.vars) < len(b.vars):
            return -1
        elif len(a.vars) > len(b.vars):
            return 1
        else:
            return a.coeff - b.coeff

    def __repr__(self):
        if self.nickname:
            return self.nickname
        else:
            return super().__repr__()


@dataclass(frozen=True)
class DegLex(TermOrdering):
    """Degree lexicographic term ordering."""

    def __call__(self, a: "Term", b: "Term") -> int:
        if a.degree() < b.degree():
            return -1
        elif a.degree() > b.degree():
            return 1
        else:
            return Lex(self.vars)(a, b)

    def __repr__(self):
        if self.nickname:
            return self.nickname
        else:
            return super().__repr__()


@dataclass(frozen=True)
class DegRevLex(TermOrdering):
    """Degree reverse lexicographic term ordering."""

    def __call__(self, a: "Term", b: "Term") -> int:
        if a.degree() < b.degree():
            return -1
        elif a.degree() > b.degree():
            return 1
        else:
            return -Lex(self.vars)(a, b)

    def __repr__(self):
        if self.nickname:
            return self.nickname
        else:
            return super().__repr__()


default_ordering = string.ascii_letters[::-1]  # ... > x > y > z
lex = Lex(default_ordering, "lex")
deglex = DegLex(default_ordering, "deglex")
degrevlex = DegRevLex(default_ordering, "degrevlex")

TERM_RE = re.compile(
    r"^\s*([+-])?\s*(\d+)?((?:\*?[a-zA-Z](?:\^[+-]?\d+)?)*)\s*([+-].*)?$"
)
VAR_RE = re.compile(r"([a-zA-Z])(?:\^([+-]?\d+))?")


@total_ordering  # all term orders are total orders
@dataclass(frozen=True)
class Term:
    """A single term in a polynomial."""

    coefficient: int
    vars: tuple[tuple[str, int], ...]
    ordering: TermOrdering

    def __post_init__(self):
        vars = defaultdict(int)
        for var, exp in self.vars:
            vars[var] += exp
            if vars[var] == 0:
                del vars[var]

        if len(vars) > 1:
            keys = sorted(
                vars.keys(), key=lambda var: Term(1, ((var, 1),), self.ordering)
            )
            object.__setattr__(
                self, "vars", tuple((key, vars[key]) for key in keys)
            )
        else:
            object.__setattr__(self, "vars", tuple(vars.items()))

    def __repr__(self):
        vars = "".join(
            var if exp == 1 else f"{var}^{exp}" for var, exp in self.vars[::-1]
        )
        if self.coefficient == 1 and vars:
            return vars
        elif self.coefficient == -1 and vars:
            return "-" + vars
        else:
            return str(self.coefficient) + vars

    @classmethod
    def parse(cls, term, ordering=lex):
        """Parse a term from a human-readable string."""
        if isinstance(term, str):
            sign, coefficient, vars, rest = TERM_RE.match(term).groups()
        else:
            sign, coefficient, vars, rest = term

        assert coefficient or vars
        if not rest:
            pass
        else:
            raise ValueError((term, sign, coefficient, vars, rest))
        assert not rest

        return Term(
            int((sign or "") + (coefficient or "1")),
            tuple(
                (m.group(1), int(m.group(2) or 1))
                for m in VAR_RE.finditer(vars)
            ),
            ordering,
        )

    def degree(self):
        """Return the degree of the term."""
        return sum(exp for _, exp in self.vars)

    def divides(self, other):
        """Check if the term would divide the other."""
        if isinstance(other, Term):
            if self.ordering != other.ordering:
                raise ValueError(
                    f"Orderings differ: {self.ordering} != {other.ordering}"
                )
            other_vars = dict(other.vars)
            return all(
                other_vars.get(var, 0) - exp >= 0 for var, exp in self.vars
            )
        elif isinstance(other, Polynomial):
            return Polynomial((self,), self.ordering).divides(other)

    def reduce(self, divisors):
        """Divide by a set of polynomials until none divide the polynomial."""
        return Polynomial((self,), self.ordering).reduce(divisors)

    def divmod(self, divisor):
        """Divide the term by another term."""
        return Polynomial((self,), self.ordering).divmod(divisor)

    def __lt__(self, other):
        if isinstance(other, Term):
            if self.ordering != other.ordering:
                raise ValueError(
                    f"Orderings differ: {self.ordering} != {other.ordering}"
                )
            return self.ordering(self, other) < 0
        elif isinstance(other, int):
            return self < Term(other, (), self.ordering)
        else:
            return NotImplemented

    def __neg__(self):
        """Negate a term."""
        return Term(-self.coefficient, self.vars, self.ordering)

    def __add__(self, other):
        """Add two terms."""
        if isinstance(other, Term):
            return Polynomial((self, other), self.ordering)
        elif isinstance(other, int):
            return self + Term(other, (), self.ordering)
        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtract two terms."""
        if isinstance(other, Term) or isinstance(other, Polynomial):
            return self + (-other)
        elif isinstance(other, int):
            return self - Term(other, (), self.ordering)
        else:
            return NotImplemented

    def __mul__(self, other):
        """Multiply two terms."""
        if isinstance(other, Term):
            return Term(
                self.coefficient * other.coefficient,
                self.vars + other.vars,
                self.ordering,
            )
        elif isinstance(other, int):
            return self * Term(other, (), self.ordering)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Term):
            if other.coefficient == 0:
                raise ZeroDivisionError(f"{self} / {other}: Division by zero")
            elif self.coefficient % other.coefficient != 0:
                raise ValueError(f"{self} / {other}: Division is not exact")
            else:
                exponents = dict(self.vars)
                for var, exp in other.vars:
                    new_exp = exponents.get(var, 0) - exp
                    if new_exp < 0:
                        raise ZeroDivisionError(
                            f"{self} / {other}: {var}^{new_exp} is not in the ring"
                        )
                    exponents[var] = new_exp
                assert (self.coefficient % other.coefficient) == 0
                return Term(
                    self.coefficient // other.coefficient,
                    tuple(exponents.items()),
                    self.ordering,
                )
        elif isinstance(other, int):
            if other == 0:
                raise ZeroDivisionError(f"{self} / {other}: Division by zero")
            elif self.coefficient % other != 0:
                raise ValueError(f"{self} / {other}: Division is not exact")
            else:
                return Term(
                    self.coefficient // other,
                    self.vars,
                    self.ordering,
                )
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, int):
            return self + other
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, int):
            return (-self) + other
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, int):
            return self * other
        else:
            return NotImplemented


@dataclass(frozen=True)
class Polynomial:
    """A polynomial in a (potentially multivariate) ring of polynomials."""

    terms: tuple[Term, ...]
    ordering: TermOrdering

    def __post_init__(self):
        coeffs = defaultdict(int)
        for term in self.terms:
            coeffs[term.vars] += term.coefficient

        terms = [
            Term(coeff, vars, self.ordering)
            for vars, coeff in coeffs.items()
            if coeff
        ]
        terms.sort()

        object.__setattr__(self, "terms", tuple(terms))

    def __str__(self):
        if self.terms:
            s = " + ".join(str(term) for term in self.terms[::-1])
            return s.replace(" + -", " - ")
        else:
            return "0"

    @classmethod
    def parse(cls, polynomial, ordering):
        """Parse a polynomial from a human-readable string representation."""

        def terms(s):
            while s:
                sign, coefficient, vars, s = TERM_RE.match(s).groups()
                yield Term.parse((sign, coefficient, vars, ""), ordering)

        return Polynomial(tuple(terms(polynomial)), ordering)

    def degree(self):
        """Return the degree of the polynomial."""
        return max(term.degree() for term in self.terms) if self.terms else 0

    def leading_term(self):
        """Return the leading term of the polynomial."""
        return self.terms[-1] if self.terms else Term(0, (), self.ordering)

    def divides(self, other):
        """Check if the polynomial would divide the other with no remainder."""
        if isinstance(other, Polynomial):
            try:
                _, remainder = other.reduce((self,))
                return not remainder.terms
            except ZeroDivisionError:
                return False
        elif isinstance(other, Term):
            return len(self.terms) == 1 and self.terms[0].divides(other)
        else:
            raise TypeError(f"Unsupported type: {type(other)}")

    def reduce(self, divisors):
        """Divide by a set of polynomials until none divide the polynomial."""

        # Multivariable Division Algorithm
        # Let f, f_1, ..., f_s in R with f_i != 0 for all i.
        # Set r := 0, h := f, u_1, ..., u_s := 0.
        # While h != 0, do
        #   If there exists i with lt(f_i) divides lt(h), then
        #       Choose i such that lt(f_i) divides lt(h).
        #       Set u_i := u_i + lt(h) / lt(f_i).
        #       Set h := h - (lt(h) / lt(f_i)) * f_i.
        #   Else
        #       Set r := r + lt(h).
        #       Set h := h - lt(h).
        # Note: At each stage of the algorithm, f = u_1f_1 + ... + u_sf_s + h + r.

        if not all(d.terms for d in divisors):
            raise ZeroDivisionError("Division by zero polynomial")

        zero = Polynomial((), self.ordering)
        quotients = {divisor: zero for divisor in divisors}
        remainder = zero

        dividend = self
        while dividend != zero:
            for divisor in divisors:
                try:
                    print(
                        f"{dividend.leading_term()} / {divisor.leading_term()}"
                    )
                    q_term = dividend.leading_term() / divisor.leading_term()
                    print(q_term)
                    print(q_term * divisor)
                    quotients[divisor] += q_term
                    dividend -= q_term * divisor
                    break
                except ZeroDivisionError:
                    pass
            else:
                remainder += dividend.leading_term()
                dividend -= dividend.leading_term()

        return quotients, remainder

    def divmod(self, divisor):
        """Divide the polynomial by another polynomial."""
        quotients, remainder = self.reduce((divisor,))
        return quotients[divisor], remainder

    def __neg__(self):
        return self * Polynomial.parse("-1", self.ordering)

    def __add__(self, other):
        """Add two polynomials."""
        if isinstance(other, Polynomial):
            if self.ordering != other.ordering:
                raise ValueError(
                    f"Orderings differ: {self.ordering} != {other.ordering}"
                )

            def terms(a, b):
                i = j = 0
                while i < len(a) and j < len(b):
                    if a[i].vars == b[j].vars:
                        coefficient = a[i].coefficient + b[j].coefficient
                        if coefficient:
                            yield Term(coefficient, a[i].vars, self.ordering)
                        i += 1
                        j += 1
                    elif a[i] < b[j]:
                        yield a[i]
                        i += 1
                    else:
                        yield b[j]
                        j += 1
                yield from a[i:]
                yield from b[j:]

            return Polynomial(
                tuple(terms(self.terms, other.terms)), self.ordering
            )
        elif isinstance(other, Term):
            return self + Polynomial((other,), self.ordering)
        elif isinstance(other, int):
            return self + Polynomial.parse(str(other), self.ordering)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Polynomial) or isinstance(other, Term):
            return self + (-other)
        elif isinstance(other, int):
            return self - Polynomial.parse(str(other), self.ordering)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            terms = defaultdict(int)
            for a in self.terms:
                for b in other.terms:
                    coeffs = defaultdict(int)
                    for var, exp in a.vars + b.vars:
                        coeffs[var] += exp

                    coefficient = a.coefficient * b.coefficient
                    terms[tuple(sorted(coeffs.items()))] += coefficient

            return Polynomial(
                tuple(
                    Term(coeff, vars, self.ordering)
                    for vars, coeff in terms.items()
                ),
                self.ordering,
            )
        elif isinstance(other, Term):
            return self * Polynomial((other,), self.ordering)
        elif isinstance(other, int):
            return self * Polynomial.parse(str(other), self.ordering)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, Polynomial):
            quotients, remainder = self.reduce((other,))
            return quotients[other]
        elif isinstance(other, Term):
            return self / Polynomial((other,), self.ordering)
        elif isinstance(other, int):
            term = Term(other, (), self.ordering)
            return self / Polynomial.parse((term,), self.ordering)
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Polynomial):
            quotients, remainder = self.reduce((other,))
            return remainder
        elif isinstance(other, Term):
            return self % Polynomial((other,), self.ordering)
        elif isinstance(other, int):
            term = Term(other, (), self.ordering)
            return self % Polynomial.parse((term,), self.ordering)
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Term) or isinstance(other, int):
            return self + other
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, Term) or isinstance(other, int):
            return (-self) + other
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Term) or isinstance(other, int):
            return self * other
        else:
            return NotImplemented
