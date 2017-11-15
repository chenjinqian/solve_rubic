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
        # for i in self.fml:
        for i in ["R"]:
            yield i

    def gen_up(self, g):
        return ("%s%s" % (i, j) for i in g
                for j in self.fml if self.check_fml_valid(i, j))

    def check_fml_valid(self, item_i, item_j, v=False):
        allow_d = {}
        for one in self.fml:
            allow_d[one] = {}
        allow_list = [
            ["R", "F", "U", "B", "D", "F'", "U'", "B'", "D'", "L", "L'"],
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
            ["R", "F", "L", "B", "R'", "F'", "L'", "B'"]]
        indep_d = {
            "R": {"L": '', "L'": ''},
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
            "D'": {"U": '', "U'": ''}}
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
            "D'": "D"}
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
                    if ((item_j == reverse_d[last_second]) or
                            (item_j == last_second)):
                        return False
        except Exception as e:
            print("%s ERROR" % (repr(e)))
            print(item_j, item_i, items)
        return True

    def gen_level(self, n):
        if n > 1:
            return self.gen_up(self.gen_level(n-1))
        return self.fml_gen1()

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

    def expand_d(self):
        pass

    def loop_from_hash(self, matrix_hash):
        """
        input:
        '6_2_3_8b_1b_4_7_12b_5_9_11_10b/14b_16c_15_20b_18c_13b_19_17c'
        output:
        (come or pull, means)
        '6.4.8b.12b.10b.9.5.1b/14b.16c.20b.17c.18c.13b'
        """
        status_d = self._status_d_from_hash(matrix_hash)
        lp = self._loop_from_status_d(status_d)
        return lp

    def _status_d_from_hash(self, matrix_hash):
        d = {}
        sybs = [i for j in matrix_hash.split("/") for i in j.split("_")]
        for num_0, syb in zip(range(20), sybs):
            d[str(int(num_0) + 1)] = syb
        return d

    def _hash_from_status_d(self, status_d):
        symbs = [status_d[str(int(i) + 1)] for i in range(20)]
        hs = "/".join(["_".join(symbs[:12]), "_".join(symbs[12:])])
        return hs

    def _loop_from_status_d(self, status_d):
        """
        need new version to speed up.
        """
        tmp_d = {}
        key_seq = ["%s" % (int(i) + 1) for i in range(20)]
        lp = []
        while (len(tmp_d) < 20):
            for key in key_seq:
                if key in tmp_d:
                    continue
                # print("#2, key %s" % (key))
                val = status_d[key]
                if (key == val):
                    # print("key=val, %s" % (key))
                    tmp_d[key] = ''
                else:
                    first_symb = val
                    # print("#4, first_symb %s" % (first_symb))
                    sub_loop = []
                    sub_loop.append(val)
                    tmp_d[key] = ''
                    if val not in status_d:
                        key = "%s" % (val[:-1])
                    else:
                        key = "%s" % (val)
                    val = status_d[key]
                    while (not val == first_symb):
                        sub_loop.append(val)
                        tmp_d[key] = ''
                        if val not in status_d:
                            key = "%s" % (val[:-1])
                        else:
                            key = "%s" % (val)
                        val = status_d[key]
                        # print("#3, val %s, key_next %s" % (val, key))
                    lp.append(sub_loop)
        lp_str = "/".join([".".join(i) for i in lp])
        return lp_str

    def _status_d_from_loop(self, lp):
        base_status_d = self._status_d_from_hash(self.hash_matrix(base))
        lp_list = [s.split(".") for s in lp.split("/")]
        for sub_loop in lp_list:
            symb_1st = ''
            place = ''
            for symb in sub_loop:
                if not symb_1st:
                    symb_1st = symb
                if place:
                    base_status_d[place] = symb
                if symb not in base_status_d:
                    symb = symb[:-1]
                place = symb
            base_status_d[place] = symb_1st
        return base_status_d

    def _status_d_operator(self, sd1, sd2):
        rlt_d = {}
        for key in sd2:
            val2 = sd2[key]
            if val2 not in sd2:
                val_num2 = val2[:-1]
                val_pos2 = val2[-1:]
            else:
                val_num2 = val2
                val_pos2 = ''
            val1 = sd1[val_num2]
            if val1 not in sd1:
                val_num1 = val1[:-1]
                val_pos1 = val1[-1:]
            else:
                val_num1 = val1
                val_pos1 = ''
            rlt_num = val_num1
            if not (val_pos1 and val_pos2):
                # ('', b), ('', c) -> b, c
                rlt_pos = (val_pos1 or val_pos2)
            elif not (val_pos1 == val_pos2):
                # (c, b) -> ''
                rlt_pos = ''
            elif int(val_num2) <= 12:
                # (b, b)_12 -> ''
                rlt_pos = ''
            elif (val_pos1 == 'b'):
                # (b, b)_20 -> 'c'
                rlt_pos = 'c'
            else:
                # (c, c)_20 -> 'b'
                rlt_pos = 'b'
            rlt_d[key] = "%s%s" % (rlt_num, rlt_pos)
        return rlt_d

    def loop_operator(self, lp1, lp2):
        sd1 = self._status_d_from_loop(lp1)
        sd2 = self._status_d_from_loop(lp2)
        sd_rlt = self._status_d_operator(sd1, sd2)
        lp_rlt = self._loop_from_status_d(sd_rlt)
        return lp_rlt

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
        rlt = "%s/%s" % ("_".join(pos_pha_list[:12]),
                         "_".join(pos_pha_list[12:]))
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

    # TODO: ipmlement one solve method of rubic cube
    # DONE: XYZ is a small space similar problem. It could help.
    # TODO: use these functinos to find out low level formula equivlent list.
    # TODO: change right multiple to left multple, for easy to write reason.
    # or just reverse the formula
    # TODO: use reduced symbol caculation, not directerly matrix multiple,
    # group operator.
    # TODO: expand_d
    # TODO: group analysis hash value, and count mct.


if __name__ == '__main__':
    main()
