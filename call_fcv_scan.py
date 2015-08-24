import sys
import os
import datetime
import tempfile
import subprocess

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('evttype', help='evttype to cross validate. 1-5.', type=int)
    parser.add_argument('colnum', help='column index cross validate. 0 or 1.', type=int)
    parser.add_argument('output_fname', help='output file name', type=str)
    parser.add_argument('--r', type=int, default=15,
                                help='Use FFT of size 2^r. ')
    parser.add_argument('--begin', type=float, default=0.1,
                                help='beginning bandwidth. ')
    parser.add_argument('--end', type=float, default=0.0001,
                                help='ending bandwidth. ')
    parser.add_argument('--step', type=float, default=0.0001,
                                help='decrement step size from start to end. ')
    parser.add_argument('--data_dir', default='data',
                                        help='Absolute path to kde input data.')
    args = parser.parse_args()

    input_fname = args.data_dir + '/evttype{0}.dither.csv'.format(args.evttype)

    subprocess.check_call(['./fcv_scan', args.output_fname, input_fname,
                           str(args.colnum), str(args.r),
                           str(args.begin), str(args.end), str(args.step)])
