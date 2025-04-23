import lhapdf
import math

Ebeam = 11
Q2  = 1.5
PhT = 0.1
xB  = 0.125
zh  = 0.4

M = 0.938
s = 2 * M * Ebeam + M**2
y = Q2 /( xB *( s - M**2 ) )
gamma = 2 * xB * M / math.sqrt(Q2)
epsilon = ( 1 - y - 1/4 * y**2 * gamma**2 )/( 1 - y + 1/2 * y**2 + 1/4 * y**2 * gamma**2 )

# one-loop single-flavor running of alphaEM
me = 0.000511
alphaEM0 = 1/137.036
def alphaEM(Q2):
    return 1/( 1/alphaEM0 - 1/(3*math.pi) * math.log(Q2/me**2) )

num_quark = 10
quark_code   = [ +2  , -2  , +1  , -1  , +3  , -3  , +4  , -4  , +5  , -5   ]
quark_charge = [ +2/3, -2/3, -1/3, +1/3, -1/3, +1/3, +2/3, -2/3, -1/3, +1/3 ]
# isospin transformation
def iso_trans(q):
    match q:
        case  2: return  1
        case  1: return  2
        case -2: return -1
        case -1: return -2
        case  _: return  q

pdf = lhapdf.mkPDF("CJ15lo/0")
ff  = lhapdf.mkPDF("dsspipLO/0")
def PhT2_avg(zh):
    return zh**2 * 0.25 + 0.2

print(
    # d sigma / dxB dy dzh dphih dPhT^2
    # ref: [Alessandro Bacchetta et al JHEP02(2007)093]
    2*math.pi *
    alphaEM(Q2)**2 /( xB * y * Q2 ) *
    y**2 /( 2*(1-epsilon) ) *
    ( 1 + gamma**2/(2*xB) ) *
    xB * sum(
        quark_charge[i]**2 *
        1/xB * 1/3 * ( 2 * pdf.xfxQ2(          quark_code[i] , xB, Q2)
                     +     pdf.xfxQ2(iso_trans(quark_code[i]), xB, Q2) ) *
        1/zh * ff.xfxQ2(quark_code[i], zh, Q2) *
        math.exp( - PhT**2 / PhT2_avg(zh) )/( math.pi * PhT2_avg(zh) )
        for i in range(0,num_quark)
    ) *
    389 * 10**3 # unit conversion to (nb/GeV^2)
)
