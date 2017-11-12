#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: chenjinqian
# email: 2012chenjinqian@gmail.com


# matrix represent  rubic state and try to solve the best path .
import numpy as np
# import scipy as sp
import cmath
import time


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
Fr = Zr * Rr * Z
F = Fr * Fr * Fr
Br = Z * Z * F * Z * Z
B = Br * Br * Br

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
        hash_d is like
        {"1_2_3_4_5b_6b_7_8_9_10_11_12/13_14_15_16_17_18_19_20":
        ["R'HR'HR2H'R'H'R'H'R2H", ...],
        ...}
        """
        self.fml = ["R", "F", "U", "L", "B", "D",
                    "R'", "F'", "U'", "L'", "B'", "D'"]
        self.d_8 = {}
        self.d_12 = {}
        self.d_20 = {}

    # def operation_in_e(self, ops):
    #     ops_revers = reversed(ops)
    #     return ops_revers

    def fml_from_count(self, cnt):
        """
        cnt is int
        """
        fml = ""
        cnt = int(cnt)
        while cnt > 0:
            low_one = cnt  % 11
            cnt = int(cnt / 12)
            fml += self.fml[low_one]
        return fml

    def fml_gen_acc(self):
        for i in self.fml:
            yield i

    def eval_fml(self, fml, left_to_right=True):
        def fml_to_operation(fml):
            operations = []
            item_a = None
            item_b = None
            for item in fml:
                item_a = item
                if item_b is None:
                    item_b = item_a
                    continue
                if item_a == "'":
                    item_b = "%s%s" % (item_b, "r")
                elif item_a == "2":
                    operations.append(item_b)
                else:
                    operations.append(item_b)
                    item_b = item_a
            operations.append(item_b)
            return operations

        def eval_operation(operations):
            """
            evalue operation;
            using left to right order.
            """
            if left_to_right:
                elements = reversed(operations)
            else:
                elements = operations
            mk = " * "
            fs = []
            for e in elements:
                if not (fs == []):
                    fs = fs + mk + e
                else:
                    fs = e
            try:
                matrix = eval(fs)
                return matrix
            except Exception as e:
                print("#2, ERROR: %s" % (repr(e)))
                return []
        return eval_operation(fml_to_operation(fml))

    def steps(self, fml):
        operations = []
        item_a = None
        item_b = None
        for item in fml:
            item_a = item
            if item_b is None:
                item_b = item_a
                continue
            if item_a == "'":
                item_b = "%s%s" % (item_b, "r")
            elif item_a == "2":
                operations.append(item_b)
            else:
                operations.append(item_b)
                item_b = item_a
        operations.append(item_b)
        return len(operations)

    def hash_matrix(self, matrix):
        """
        f2: R'HR'HR2H'R'H'R'H'R2H
        f2 hash: 1_2_3_4_5b_6b_7_8_9_10_11_12/13_14_15_16_17_18_19_20
        """
        # abc_str = "abcdefghijklmnopqrst"
        # num_abc_d = {}
        # for num, abc in zip(range(20), abc_str):
        #     num_abc_d["%s" % (num)] = abc

        def position(matrix):
            """
            should transform before get shape.
            """
            p_019 = list(matrix.T.nonzero()[1])
            p_120 = [i + 1 for i in p_019]
            return p_120

        def phase(matrix):
            m201 = c_20_1 * matrix
            clx_angle = list(m201.A1)
            # print("clx_angle %s" % (clx_angle))
            phase_pi = [cmath.phase(i) for i in clx_angle]
            digt = []
            for i in phase_pi:
                if i > 0.1:
                    digt.append('b')
                elif i < -0.1:
                    digt.append('c')
                else:
                    digt.append('a')
            return digt

        pos_pha_list = []
        for pos, pha in zip(position(matrix), phase(matrix)):
            if "a" == "%s" % (pha):
                pos_pha_list.append("%s" % (pos))
            else:
                pos_pha_list.append("%s%s" % (pos, pha))
        rlt = "%s/%s" % ("_".join(pos_pha_list[:12]), "_".join(pos_pha_list[12:]))
        return rlt


def main():
    f1 = "R'U'2RUR'URURU2R'U'RU'R'"
    f2 = "R'HR'HR2H'R'H'R'H'R2H"
    f3 = "HR2H'R'H'R'H'R2HR'HR'"
    s = RubicMatrix()
    print("f1: %s, %s steps." % (f1, s.steps(f1)))
    print("f2: %s, %s steps." % (f2, s.steps(f2)))
    print("f3: %s, %s steps." % (f3, s.steps(f3)))
    ta = time.time()
    print("f1 hash: %s" % (s.hash_matrix(s.eval_fml(f1))))
    print("f2 hash: %s" % (s.hash_matrix(s.eval_fml(f2))))
    print("f3 hash: %s" % (s.hash_matrix(s.eval_fml(f3))))
    tb = time.time()
    print("spend %s s" % (float(tb - ta)))
    # TODO: ipmlement one solve method of rubic cube
    # DONE: XYZ is a small space similar problem. It could help.
    # TODO: use these functinos to find out low level formula equivlent list.
    # TODO: change right multiple to left multple, for easy to write reason.
    # or just reverse the formula


if __name__ == '__main__':
    main()
