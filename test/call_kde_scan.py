import sys
import os
import datetime
import tempfile
import subprocess

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputf', type=str,
                        help='Sample points to form the kde. ')
    parser.add_argument('h1', type=float,
                        help='Bandwidth in first dimension. ')
    parser.add_argument('h2', type=float,
                        help='Bandwidth in second dimension. ')
    parser.add_argument('--outputf_2d', type=str, default='prodkde2d_scan2d.dat',
                        help='File name to save 2d scan. ')
    parser.add_argument('--outputf_1d', type=str, default='prodkde2d_scan1d.dat',
                        help='File name to save 1d scan. ')
    parser.add_argument('--b1', type=float, default=-1.5,
                        help='Starting scan point for the first dimension. ')
    parser.add_argument('--e1', type=float, default=2.0,
                        help='End point bound for the first dimension. ')
    parser.add_argument('--s1', type=float, default=0.05,
                        help='Step size in the first dimension. ')
    parser.add_argument('--b2', type=float, default=-0.5,
                        help='Starting scan point for the second dimension. ')
    parser.add_argument('--e2', type=float, default=2.0,
                        help='End point bound for the second dimension. ')
    parser.add_argument('--s2', type=float, default=0.01,
                        help='Step size in the second dimension. ')
    args = parser.parse_args()


    subprocess.check_call(['./prodkde2d_scan',
                           args.outputf_2d, args.outputf_1d, args.inputf,
                           str(args.h1), str(args.h2),
                           str(args.b1), str(args.e1), str(args.s1),
                           str(args.b2), str(args.e2), str(args.s2)])
