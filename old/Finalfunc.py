from Phase1 import phase1
from Phase2 import phase2
from Phase3 import phase3

def finalfunc(h, s, bc):
    p = [9.81, 12, 4000, 100, 8, h, s, bc, 6.5]
    S1 = phase1(p)
    S2 = phase2(p, S1)
    S3 = phase3(p, S2)
    return S1, S2, S3
    
a, b, c = finalfunc(6, 6, 1.5)