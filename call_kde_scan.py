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
# naive cv
#opt_h1 = [ 2.29e-3, 2.10e-3, 1.24e-3, 1.42e-3, 1.65e-4 ]
#opt_h2 = [ 2.16e-2, 1.64e-2, 1.91e-2, 1.63e-2, 1.44e-2 ]
# fft cv
opt_h1 = [ 8.76e-3, 7.33e-3, 2.81e-3, 2.92e-3, 1.49e-5 ]
opt_h2 = [ 7.84e-2, 8.04e-2, 8.05e-2, 5.41e-2, 3.85e-2 ]

for i, (h1, h2) in enumerate(zip(opt_h1, opt_h2)):
    sys.stdout.write('Scanning evttype{0}...\n'.format(i+1))
    sys.stdout.flush()
    out2d = outdir + '/evttype{0}.2d.fcv.opt.txt'.format(i+1)
    out1d = outdir + '/evttype{0}.1d.fcv.opt.txt'.format(i+1)
    in_data = datadir + '/evttype{0}.dither.csv'.format(i+1)

    subprocess.check_call(['./kde_scan', out2d, out1d, in_data,
                           str(h1), str(h2),
                           str(b1), str(e1), str(s1),
                           str(b2), str(e2), str(s2)])
