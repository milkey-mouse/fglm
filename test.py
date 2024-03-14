#!/usr/bin/env python3
import unittest
from polynomials import *
from groebner import *


class TestTerm(unittest.TestCase):
    def test_term_init(self):
        t1 = Term.parse("2x^2y", lex)
        self.assertEqual(str(t1), "2x^2y")

    def test_term_parse(self):
        t1 = Term.parse("2x^2y", lex)
        self.assertEqual(str(t1), "2x^2y")

    def test_term_degree(self):
        t1 = Term.parse("2x^2y", lex)
        self.assertEqual(t1.degree(), 3)

    def test_divides(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("xy", lex)
        self.assertTrue(t2.divides(t1))
        self.assertFalse(t1.divides(t2))

    def test_reduce(self):
        t1 = Term.parse("2x^2y", lex)
        p1 = Polynomial.parse("x + y", lex)
        print(t1.reduce((p1,)))
        quotients, remainder = t1.reduce((p1,))
        self.assertEqual(str(quotients[p1]), "2x^2y")

    def test_term_lt(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("3xy^2", lex)
        self.assertLess(t1, t2)

    def test_term_neg(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = -t1
        self.assertEqual(str(t2), "-2x^2y")

    def test_term_add(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("3xy^2", lex)
        p = t1 + t2
        self.assertEqual(str(p), "3xy^2 + 2x^2y")

    def test_term_sub(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("x^2y", lex)
        p = t1 - t2
        self.assertEqual(str(p), "x^2y")

    def test_term_mul(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("3xy^2", lex)
        t3 = t1 * t2
        self.assertEqual(str(t3), "6x^3y^3")

    def test_term_truediv(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("xy", lex)
        t3 = t1 / t2
        self.assertEqual(str(t3), "2x")

    def test_term_radd(self):
        t1 = Term.parse("2x^2y", lex)
        p = 3 + t1
        self.assertEqual(str(p), "2x^2y + 3")

    def test_term_rsub(self):
        t1 = Term.parse("2x^2y", lex)
        p = 3 - t1
        self.assertEqual(str(p), "-2x^2y + 3")

    def test_term_rmul(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = 3 * t1
        self.assertEqual(str(t2), "6x^2y")


class TestPolynomial(unittest.TestCase):
    def test_polynomial_init(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        self.assertEqual(str(p1), "2x^2y + 3y")

    def test_polynomial_parse(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        self.assertEqual(str(p1), "2x^2y + 3y")

    def test_polynomial_degree(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        self.assertEqual(p1.degree(), 3)

    def test_polynomial_leading_term(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        self.assertEqual(str(p1.leading_term()), "2x^2y")

    def test_polynomial_reduce(self):
        p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
        p2 = Polynomial.parse("x + y", lex)
        quotients, remainder = p1.reduce((p2,))
        self.assertEqual(str(quotients[p2]), "x + y")
        self.assertEqual(str(remainder), "0")

    def test_reduce_single_variable(self):
        f = Polynomial.parse("x^3 + x^2 + x + 1", lex)
        f1 = Polynomial.parse("x^2 + 1", lex)
        f2 = Polynomial.parse("x + 1", lex)

        quotients, remainder = f.reduce((f1, f2))
        self.assertEqual(str(quotients[f1]), "x + 1")
        self.assertEqual(str(quotients[f2]), "x^2 - 1")
        self.assertEqual(str(remainder), "0")

    def test_reduce_multiple_variables(self):
        f = Polynomial.parse("x^2*y + x*y^2 + y", lex)
        f1 = Polynomial.parse("x*y - 1", lex)
        f2 = Polynomial.parse("y^2 - 1", lex)

        quotients, remainder = f.reduce((f1, f2))
        self.assertEqual(str(quotients[f1]), "x + y")
        self.assertEqual(str(quotients[f2]), "x")
        self.assertEqual(str(remainder), "x + y + 1")

    def test_reduce_different_orderings(self):
        f = Polynomial.parse("x^2*y + x*y^2 + y", lex)
        f1 = Polynomial.parse("x*y - 1", deglex)
        f2 = Polynomial.parse("y^2 - 1", degrevlex)

        with self.assertRaises(ValueError):
            f.reduce((f1, f2))

    def test_reduce_zero_divisor(self):
        f = Polynomial.parse("x^2 + x + 1", lex)
        f1 = Polynomial.parse("0", lex)

        with self.assertRaises(ZeroDivisionError):
            f.reduce((f1,))

    def test_reduce_no_division(self):
        f = Polynomial.parse("x^2 + x + 1", lex)
        f1 = Polynomial.parse("y^2 + y + 1", lex)

        quotients, remainder = f.reduce((f1,))
        self.assertEqual(str(quotients[f1]), "0")
        self.assertEqual(str(remainder), "x^2 + x + 1")

    def test_reduce_multiple_divisors(self):
        f = Polynomial.parse("x^3*y^2 + x^2*y^3 + x*y^4 + y^5", lex)
        f1 = Polynomial.parse("x*y^2 - 1", lex)
        f2 = Polynomial.parse("y^3 - 1", lex)

        quotients, remainder = f.reduce((f1, f2))
        self.assertEqual(str(quotients[f1]), "x^2 + x*y + y^2")
        self.assertEqual(str(quotients[f2]), "x*y + y^2")
        self.assertEqual(str(remainder), "x^2 + x*y + y^2 + 1")

    def test_reduce_large_polynomial(self):
        f = Polynomial.parse(
            "a^5*b^4*c^3*d^2*e + a^4*b^3*c^2*d*e^2 + a^3*b^2*c*d^3*e^3 + a^2*b*c^4*d^4*e^4 + a*b^5*c^5*d^5*e^5",
            lex,
        )
        f1 = Polynomial.parse("a*b*c*d*e - 1", lex)
        f2 = Polynomial.parse("a^2*b^2*c^2*d^2*e^2 - 1", lex)

        quotients, remainder = f.reduce((f1, f2))
        self.assertEqual(
            str(quotients[f1]),
            "a^4*b^3*c^2*d + a^3*b^2*c*d^2*e + a^2*b*c^3*d^3*e^2 + a*b^4*c^4*d^4*e^3 + b^5*c^5*d^5*e^4",
        )
        self.assertEqual(
            str(quotients[f2]),
            "a^3*b^2*c + a^2*b*c^2*d + a*b^3*c^3*d^2*e + b^4*c^4*d^3*e^2",
        )
        self.assertEqual(
            str(remainder),
            "a^4*b^3*c^2*d + a^3*b^2*c*d^2*e + a^2*b*c^3*d^3*e^2 + a*b^4*c^4*d^4*e^3 + b^5*c^5*d^5*e^4 + 1",
        )

    def test_divides(self):
        p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
        p2 = Polynomial.parse("x + y", lex)
        self.assertTrue(p2.divides(p1))
        self.assertFalse(p1.divides(p2))

        t1 = Term.parse("x^2", lex)
        self.assertTrue(t1.divides(p1))
        self.assertFalse(t1.divides(p2))

    def test_polynomial_divide_by(self):
        p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
        p2 = Polynomial.parse("x + y", lex)
        quotient, remainder = p1.divmod(p2)
        self.assertEqual(str(quotient), "x + y")
        self.assertEqual(str(remainder), "0")

    def test_polynomial_neg(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = -p1
        self.assertEqual(str(p2), "-2x^2y - 3y")

    def test_polynomial_add(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = Polynomial.parse("x^2y + 2y", lex)
        p3 = p1 + p2
        self.assertEqual(str(p3), "3x^2y + 5y")

    def test_polynomial_sub(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = Polynomial.parse("x^2y + 2y", lex)
        p3 = p1 - p2
        self.assertEqual(str(p3), "x^2y + y")

    def test_polynomial_mul(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = Polynomial.parse("x + 2y", lex)
        p3 = p1 * p2
        self.assertEqual(str(p3), "4x^2y^2 + 6y^2 + 2x^3y + 3xy")

    def test_polynomial_truediv(self):
        p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
        p2 = Polynomial.parse("x + y", lex)
        quotient = p1 / p2
        self.assertEqual(str(quotient), "x + y")

    def test_polynomial_mod(self):
        p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
        p2 = Polynomial.parse("x + y", lex)
        remainder = p1 % p2
        self.assertEqual(str(remainder), "0")

    def test_polynomial_radd(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = 2 + p1
        self.assertEqual(str(p2), "2x^2y + 3y + 2")

    def test_polynomial_rsub(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = 2 - p1
        self.assertEqual(str(p2), "-2x^2y - 3y + 2")

    def test_polynomial_rmul(self):
        p1 = Polynomial.parse("2x^2y + 3y", lex)
        p2 = 2 * p1
        self.assertEqual(str(p2), "4x^2y + 6y")


class TestTermOrdering(unittest.TestCase):
    def test_lex(self):
        t1 = Term.parse("2x^2y", lex)
        t2 = Term.parse("3xy^2", lex)
        self.assertLess(t1, t2)

    def test_deglex(self):
        t1 = Term.parse("2x^2y", deglex)
        t2 = Term.parse("3xy^2", deglex)
        self.assertLess(t1, t2)

    def test_degrevlex(self):
        t1 = Term.parse("2x^2y", degrevlex)
        t2 = Term.parse("3xy^2", degrevlex)
        self.assertLess(t1, t2)


class TestGroebner(unittest.TestCase):
    # def test_lcm(self):
    #    t1 = Term.parse("2x^2y", lex)
    #    t2 = Term.parse("3xy^2", lex)
    #    t3 = lcm(t1, t2)
    #    self.assertEqual(str(t3), "6x^2y^2")
    #
    # def test_S(self):
    #    p1 = Polynomial.parse("x^2 + 2xy + y^2", lex)
    #    p2 = Polynomial.parse("x + y", lex)
    #    s = S(p1, p2)
    #    self.assertEqual(str(s), "xy + y^2 - x - y")
    #
    # def test_buchberger(self):
    #    p1 = Polynomial.parse("x^2 - y^2", lex)
    #    p2 = Polynomial.parse("x^2 + y^2 - 1", lex)
    #    G = buchberger([p1, p2])
    #    self.assertEqual(len(G), 3)
    #    self.assertEqual(str(G[0]), "x")
    #    self.assertEqual(str(G[1]), "y^2 - 1/2")
    #    self.assertEqual(str(G[2]), "y")

    def test_quotient_monomials(self):
        p1 = Polynomial.parse("x^2 - y^2", lex)
        p2 = Polynomial.parse("x^2 + y^2 - 1", lex)
        G = buchberger([p1, p2])
        monomials = list(quotient_monomials(G))
        self.assertEqual(len(monomials), 4)
        self.assertEqual(str(monomials[0]), "1")
        self.assertEqual(str(monomials[1]), "y")
        self.assertEqual(str(monomials[2]), "x")
        self.assertEqual(str(monomials[3]), "xy")

    # def test_fglm(self):
    #    p1 = Polynomial.parse("x^2 - y^2", lex)
    #    p2 = Polynomial.parse("x^2 + y^2 - 1", lex)
    #    G = buchberger([p1, p2])
    #    B = fglm(G, deglex)
    #    self.assertEqual(len(B), 3)
    #    self.assertEqual(str(B[0]), "x")
    #    self.assertEqual(str(B[1]), "y^2 - 1/2")
    #    self.assertEqual(str(B[2]), "y")


if __name__ == "__main__":
    unittest.main()
