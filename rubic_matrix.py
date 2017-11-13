#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: chenjinqian
# email: 2012chenjinqian@gmail.com


# matrix represent  rubic state and try to solve the best path .
import numpy as np
# import scipy as sp
import cmath
import time
# import itertools


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
Xr = Z * Yr * Zr
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
B = Z * Z * F * Z * Z
Br = B * B * B

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
        self.fml = [
            "R", "F", "U", "L", "B", "D",
            "R'", "F'", "U'", "L'", "B'", "D'"
        ]
        self.d_8 = {}
        self.d_12 = {}
        self.d_20 = {}

    def fml_gen1(self):
        for i in ["R"]:
        # for i in self.fml:
            yield i

    def gen_up(self, g):
        # g2 = ("%s%s"%(i, j)  for i in s.fml_gen1() for j in s.fml if not ((j+"'" == i[-1]) or (i[-1]+"'" == j)))
        return ("%s%s"%(i, j) for i in g for j in self.fml  if self.check_fml_valid(i, j))

    def check_fml_valid(self, item_i, item_j, v=False):
        # if ((item_j+"'" == item_i[-1]) or (item_i[-1]+"'" == item_j)):
        #     return False
        allow_d = {}
        for one in self.fml:
            allow_d[one] = {}
        allow_list = [["R", "F", "U", "B", "D", "F'", "U'", "B'", "D'", "L", "L'"],
                      ["F", "R", "U", "L", "D", "R'", "U'", "L'", "D'", "B", "B'"],
                      ["U", "R", "F", "L", "B", "R'", "F'", "L'", "B'", "D", "D'"],
                      ["L", "F", "U", "B", "D", "R'", "F'", "U'", "B'", "D'"],
                      ["B", "R", "U", "L", "D", "R'", "F'", "U'", "L'", "D'"],
                      ["D", "R", "F", "L", "B", "R'", "F'", "U'", "L'", "B'"],
                      ["F", "U", "B", "D", "F'", "U'", "L'", "B'", "D'"],
                      ["R", "U", "L", "D", "R'", "U'", "L'", "B'", "D'"],
                      ["R", "F", "L", "B", "R'", "F'", "L'", "B'", "D'"],
                      ["F", "U", "B", "D", "F'", "U'", "B'", "D'"],
                      ["R", "U", "L", "D", "R'", "U'", "L'", "D'"],
                      ["R", "F", "L", "B", "R'", "F'", "L'", "B'"],]
        indep_d = {"R": {"L": '', "L'": ''},
                   "F": {"B": '', "B'": ''},
                   "B": {"F": '', "F'": ''},
                   "U": {"D": '', "D'": ''},
                   "L": {"R": '', "R'": ''},
                   "D": {"U": '', "U'": ''},
                   "R'": {"L": '', "L'": ''},
                   "F'": {"B": '', "B'": ''},
                   "B'": {"F": '', "F'": ''},
                   "U'": {"D": '', "D'": ''},
                   "L'": {"R": '', "R'": ''},
                   "D'": {"U": '', "U'": ''}
        }
        reverse_d = {
            "R": "R'",
            "F": "F'",
            "B": "B'",
            "U": "U'",
            "L": "L'",
            "D": "D'",
            "R'": "R",
            "F'": "F",
            "B'": "B",
            "U'": "U",
            "L'": "L",
            "D'": "D",
        }
        for first, possible_list in zip(self.fml, allow_list):
            for possible_one in possible_list:
                allow_d[first][possible_one] = ''
        items = []
        a = ''
        b = ''
        try:
            for one in item_i:
                b = one
                if "'" == b:
                    a = "%s%s" % (a, b)
                    items.append(a)
                else:
                    if a:
                        items.append(a)
                    a = b
            if not ("'" == b):
                items.append(b)
            if v:
                print("items %s" % (items))
            last_one = items[-1]
            if item_j not in allow_d[last_one]:
                return False
            if len(items) > 1:
                last_second = items[-2]
                if v:
                    print("item_i: %s, item_j: %s, last: %s, second: %s " % (
                        item_i, item_j, last_one, last_second
                    ))
                if last_one == last_second:
                    return False
                if item_j in indep_d[last_one]:
                    if (item_j == reverse_d[last_second]) or (item_j == last_second):
                        return False
        except Exception as e:
            print("%s ERROR" % (repr(e)))
            print(item_j, item_i, items)
        return True


    def gen_level(self, n):
        if n > 1:
            return self.gen_up(self.gen_level(n-1))
        return self.fml_gen1()

    # def fml_from_count(self, cnt):
    #     """
    #     cnt is int
    #     """
    #     fml = ""
    #     cnt = int(cnt)
    #     while cnt > 0:
    #         low_one = cnt  % 11
    #         cnt = int(cnt / 12)
    #         fml += self.fml[low_one]
    #     return fml

    # def fml_gen_acc(self, fml, cnt):
    #     """
    #     really hard to write in python, without macro.
    #     """
    #     print(fml, cnt)
    #     if cnt > 0:
    #         for i in self.fml:
    #             self.fml_gen_acc("%s%s" % (fml, i), cnt -1)
    #     for i in self.fml:
    #         yield "%s%s" % (fml, i)

    # def fml_gen4(self):
    #     for i in self.fml:
    #         for j in self.fml:
    #             if ((j+"'" == i) or (i+"'" == j)):
    #                 continue
    #             for k in self.fml:
    #                 if ((j+"'" == k) or (k+"'" == j)):
    #                     continue
    #                 for m in self.fml:
    #                     if ((j+"'" == k) or (k+"'" == j)):
    #                         continue
    #                     yield "%s%s%s%s" % (i, j, k, m)

    def fill_up_d(self, n):
        gen_n = self.gen_level(n)
        cnt = 0
        for fml in gen_n:
            try:
                hash_20 = self.hash_matrix(self.eval_fml(fml))
                hash_12, hash_8 = hash_20.split("/")
                # check_point = hash_8.split("_")[-1]
                # if not "20" == check_point:
                #     continue
                # 1 / 24 situation need to be considered.
                if hash_20 in self.d_20:
                    self.d_20[hash_20] = self.d_20[hash_20] + [fml]
                else:
                    self.d_20[hash_20] = [fml]
                if hash_12 in self.d_12:
                    self.d_12[hash_12] = self.d_12[hash_12] + [fml]
                else:
                    self.d_12[hash_12] = [fml]
                if hash_8 in self.d_8:
                    self.d_8[hash_8] = self.d_8[hash_8] + [fml]
                else:
                    self.d_8[hash_8] = [fml]
                cnt += 1
            except Exception as e:
                print(repr(e))
                print(hash_20, fml)
        return cnt

    def check_cofflict(self, d):
        for key in d:
            val = d[key]
            if len(val) > 1:
                yield (key, val)

    def eval_fml(self, fml):
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
    # TODO: use reduced symbol caculation, not directerly matrix multiple,
    # group operator.


if __name__ == '__main__':
    main()
