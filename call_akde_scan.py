import sys
import os
import datetime
import tempfile
import subprocess

b1, e1, s1 = 0.0, 0.1, 0.001
b2, e2, s2 = 0.0, 1.0, 0.01

outdir = 'visual/contour'
datadir = 'data'

# optimal cv bandwidths
# fft cv
opt_h1 = [ 1.97e-3, 1.75e-3, 3.92e-4, 5.75e-4, 1.09e-4 ]
opt_h2 = [ 1.95e-2, 1.23e-2, 6.46e-3, 7.66e-3, 4.33e-3 ]
a1, a2 = 0.5, 0.5

for i, (h1, h2) in enumerate(zip(opt_h1, opt_h2)):
    sys.stdout.write('Scanning evttype{0}...\n'.format(i+1))
    sys.stdout.flush()
    out2d = outdir + '/evttype{0}.2d.ada.txt'.format(i+1)
    out1d = outdir + '/evttype{0}.1d.ada.txt'.format(i+1)
    in_data = datadir + '/evttype{0}.dither.csv'.format(i+1)

    subprocess.check_call(['./akde_scan', out2d, out1d, in_data,
                           str(h1), str(h2),
                           str(a1), str(a2), '20',
                           str(b1), str(e1), str(s1),
                           str(b2), str(e2), str(s2)])
