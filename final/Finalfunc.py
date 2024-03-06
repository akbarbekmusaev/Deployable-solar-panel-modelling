from Phase1 import phase1
from Phase2 import phase2
from Phase3 import phase3
from Trebuchet_Functions import trebuchet_demo

def finalfunc(h, s, bc):
    p = [9.81, 12, 5000, 200, 12, h, s, bc, 12 - bc]
    S1 = phase1(p)
    S2 = phase2(p, S1)
    S3 = phase3(p, S2)
    return S1, S2, S3
