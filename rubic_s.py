#!/usr/bin/env python


HEX_LST = "1,2,3,4,5,6,7,8,9,a,b,c".split(",")
HEX_D = {"1":1,"2":2,"3":3,"4":4,
         "5":5,"6":6,"7":7,"8":8,
         "9":9,"a":10,"b":11,"c":12}
RPD_REV = {"_":"_","+":"-","-":"+",":":":"}
RPD = {"_:":":",":_":":","::":"_",
           "_+":"+","_-":"-","__":"_",
           "++":"-","+-":"_","+_":"+",
           "-+":"_","--":"+","-_":"-"}

def operator():
    d = {}
    d["E"] = "1_2_3_4_5_6_7_8_|1_2_3_4_5_6_7_8_9_a_b_c_"
    d["Z"] = "4_1_2_3_8_5_6_7_|4_1_2_3_8_5_6_7_c_9_a_b_"
    d["R"] = "5-1+3_4_6+2-7_8_|1_5:3_4_a:2_7_8_9_6_b_c_"
    d["r"] = ""
    d["F"] = ""
    d["f"] = ""
    d["B"] = ""
    d["b"] = ""
    d["U"] = ""
    d["u"] = ""
    d["D"] = ""
    d["d"] = ""
    d["L"] = ""
    d["l"] = ""
    return d

def o(s, m):
    # s "1_2_3_4_5_6_7_8_|1_2_3_4_5_6_7_8_9_a_b_c_"
    # m "5-1+3_4_6+2-7_8_|1_5:3_4_a:2_7_8_9_6_b_c_"
    rpd = RPD
    hex_d =  HEX_D
    s1, s2 = s.split("|")
    m1, m2 = m.split("|")
    len_s1m1 = min(int(len(s1)/2), int(len(m1)/2))
    len_s2m2 = min(int(len(s2)/2), int(len(m2)/2))
    rlt = "|".join([
        "".join(["%s%s"%(s1[hex_d[m1[i*2]] * 2 - 2],
                         rpd["%s%s"%(
                             s1[hex_d[m1[i*2]] * 2 - 1],
                             m1[i*2 + 1])])
                 for i in range(len_s1m1)]),
        "".join(["%s%s"%(s2[hex_d[m2[i*2]] * 2 - 2],
                         rpd["%s%s"%(
                             s2[hex_d[m2[i*2]] * 2 - 1],
                             m2[i*2 + 1])])
                 for i in range(len_s2m2)])])
    return rlt



def r(s):
    hex_lst = HEX_LST
    rpd_rev = RPD_REV
    s1, s2 = s.split("|")
    len_s1 = int(len(s1)/2)
    len_s2 = int(len(s2)/2)
    s1_mx = [[s1[i*2:i*2+2],
              "%s%s"%(hex_lst[i],
                      rpd_rev["%s"%(s1[i*2+1])])]
             for i in range(len_s1)]
    s2_mx = [[s2[i*2:i*2+2],
              "%s%s"%(hex_lst[i],
                      rpd_rev["%s"%(s2[i*2+1])])]
             for i in range(len_s2)]
    rlt = "|".join(["".join([i[1] for i in sorted(s1_mx)]),
                    "".join([i[1] for i in sorted(s2_mx)])])
    return rlt


def search(cnt=4):
    def _s(a, b, c):
        pass

    rlt = None
    return rlt



def main():
    pass

if __name__ == '__main__':
    main()
