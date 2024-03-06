from Trebuchet_Functions import trebuchet_demo
from Trebuchet_Functions import animate_trebuchet


p = [9.81, 12, 4000, 100, 8, 4, 6.4, 3.3, 4.7]

g       = p[0]
m_B     = p[1]
M_CW    = p[2]
M_P     = p[3]
L_B     = p[4]
H       = p[5]
L_S     = p[6]
L_BC    = p[7]
L_BP    = p[8]


S1, S2, S3 = trebuchet_demo(4, 6.4, 3.3)
ani = animate_trebuchet(S1.t, S1.y, S2.t, S2.y, S3.t, S3.y, p)

