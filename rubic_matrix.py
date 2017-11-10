#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: chenjinqian
# email: 2012chenjinqian@gmail.com


# matrix represent  rubic state and try to solve the best path .
import numpy as np
# import scipy as sp
import cmath


def t(n):
    return 1 * cmath.exp(1j * ((2/3.0) * cmath.pi * n))
# TODO: use decimal here.
# TODO: draw text picture to show space position.
# position is like upside to downside, side to coner,
# clock rotation, near block first.
#       position
# block
a = t(1)
b = t(2)
c = 0
d = -1
BASE = [
    [1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1],
    ]
# up reverse basic-representatio
Ur_raw = [
    [c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1],
    ]
Rr_raw = [
    [1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, d, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, d, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, a, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, b],
    [c, c, c, c, c, c, c, c, c, c, c, c, b, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, a, c, c, c],
    ]
Zr_raw = [
    [c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, 1, c, c, c],
    ]
Yr_raw = [
    [c, c, c, c, d, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, d, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, d, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, d, c, c, c, c, c, c, c, c, c, c, c],
    [1, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, d, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, 1, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, d, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, d, c, c, c, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, d, c, c, c, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, b, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, a, c, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, b, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, a],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, a, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, b, c, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, a, c, c, c, c, c],
    [c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, c, b, c],
    ]
C_20_1 = [1 for i in range(20)]
c_20_1 = np.matrix(C_20_1)
base = np.matrix(BASE)
# basic operation define.
Y = np.matrix(Yr_raw)
Yr = Y * Y * Y
Zr = np.matrix(Zr_raw)
Z = Zr * Zr * Zr
Xr = Z * Y * Zr
X = Xr * Xr * Xr
Ur = np.matrix(Ur_raw)
U = Ur * Ur * Ur
Dr = Y * Y * Ur * Y * Y
D = Dr * Dr * Dr
Rr = np.matrix(Rr_raw)
R = Rr * Rr * Rr
Lr = Z * Z * R * Z * Z
L = Lr * Lr * Lr
Br = Zr * R * Z
B = Br * Br * Br
Fr = Z * Z * B * Z * Z
F = Fr * Fr * Fr

# two layer operation, will change center orientation, not recommanded.
u = Z * D
ur = u*u*u
d = Zr * U
dr = d*d*d
r = X * L
rr = r*r*r
l = Xr * R
lr = l * l * l
f = Y * B
fr = f*f*f
b = Yr * F
br = b*b*b
# horizon, as z-ori middle layer rotation.
H = u * Ur
Hr = ur * U
# verticall, y-ori
V = f * Fr
Vr = fr * F
# slice,  x-ori
S = r * Rr
Sr = rr * R
s3 = R * u * R * U * R * U * R * u * r * u * r * R
mys3 = R * U * U * r * u * R * u * r
# too suppress too much zero after point
np.set_printoptions(suppress=True)
f22 = R*u*u*R*R*u*R*R*u*R*R*u*u*R
# colomn will be slot, row as element.
# Ur/Rr as CounterClock roatation.
# r or u as two portation rotation.
# alias
# if operation number distrane is times of 3, then
# this two operation interactive indifferent.
# times of 6, inverse operation.
e1 = R
e2 = F
e3 = U
e4 = L
e5 = B
e6 = D
e7 = Rr
e8 = Fr
e9 = Ur
e10 = Lr
e11 = Br
e12 = Dr


class RubicMatrix(object):
    def __init__(self):
        """
        formula is like   R'U'2
        operation is like  RrUU
        """
        decompose_d = {}
        self.dc = decompose_d

    # never write R2 as R'2, since it can not recongnize this form yet.
    def fml_to_operation(self, fml):
        operations = []
        a = None
        b = None
        for item in fml:
            a = item
            if b is None:
                b = a
                continue
            if a == "'":
                b = "%s%s" % (b, "r")
            elif a == "2":
                operations.append(b)
            else:
                operations.append(b)
                b = a
        operations.append(b)
        return operations

    def operation_in_e(self, ops):
        ops_revers = reversed(ops)
        return ops_revers

    def ep(self, operations):
        """
        evalue operation
        """
        elements = operations
        mk = '*'
        fs = []
        for e in elements:
            if not (fs == []):
                fs = fs + mk + e
            else:
                fs = e
        try:
            matrix = eval(fs)
            return matrix
        except:
            return []

    def ef(self, fml):
        return self.ep(self.fml_to_operation(fml))

    def position(self, matrix):
        return list(matrix.nonzero()[1])

    def oritation(self, matrix):
        return (c_20_1 * matrix)




f2 = "R'HR'HR2H'R'H'R'H'R2H"
# f3  = 'RR'
# f4 = "USSUSSUSSUSS"
# print(ef(f2).nonzero())
# print(ef(f2))
# print(ef(f3).nonzero())
# print(ef(f4).nonzero())
# print(ef("R'").nonzero())
s = RubicMatrix()
print("f2 position: %s" % (s.position(s.ef(f2))))
print("f2 oritation: %s" % (s.oritation(s.ef(f2))))
print("e1 position: %s" % (s.position(e1)))
print("f2 oritation: %s" % (s.oritation(e1)))
# TODO: simple_show, deconpose
