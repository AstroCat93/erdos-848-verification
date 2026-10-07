import sys, math
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analytic'))
import rigorous as R
from fractions import Fraction as Fr
# gain for one unsplit q-class (q=101, x odd, b in O odd) as a function of N
def gain(N, q=101, div=2):
    Lq = q*q*div
    st = {q: ('never',), 2: ('sub', Fr(1,2)), 5: ('sub', Fr(1,25))}
    one = R.class_bound(N, Lq, st, 13, 1000)
    return Fr(-(-N // Lq)) - one, one
N0 = 42_000_000
prev = None; drops = 0; ex = None
for N in range(N0, N0 + 200000, 1):
    g, one = gain(N)
    if prev is not None and g < prev[0]:
        drops += 1
        if ex is None: ex = (N - 1, float(prev[0]), N, float(g))
    prev = (g, one)
print('gain(q=101) decreases between consecutive N:', drops, 'times in 2e5 steps; first example', ex)
# Robin constant: omega(N9)*lnln/ln for N9 = primorial(23)
N9 = 2*3*5*7*11*13*17*19*23
print('Robin c at primorial(23):', 9*math.log(math.log(N9))/math.log(N9))
