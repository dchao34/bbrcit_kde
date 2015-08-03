import sys
import os
import datetime
import tempfile
import subprocess

def write_condor_submit_file(
    f, dirname, data_dir, 
    evttype, colnum, begin, end, step):

  data_dir = os.getcwd() + '/' + data_dir
  dirname = os.getcwd() + '/' + dirname
  result_fname = dirname + '/evttype{0}.col{1}.txt'.format(evttype, colnum)
  sample_fname = data_dir + '/evttype{0}.cv.csv'.format(evttype)
  out_fname = '{0}/job.out'.format(dirname)
  err_fname = '{0}/job.err'.format(dirname)
  log_fname = '{0}/job.log'.format(dirname)

  f.write('executable = cv_scan\n')
  f.write('universe = vanilla\n')
  f.write('arguments = {0} {1} {2} {3} {4} {5}\n'.format(result_fname, sample_fname, colnum, 
                                                         begin, end, step))
  f.write('output = {0}\n'.format(out_fname))
  f.write('error = {0}\n'.format(err_fname))
  f.write('log = {0}\n'.format(log_fname))
  f.write('accounting_group = group_babar\n')
  f.write('accounting_group_user = dchao\n')
  f.write('getenv = True\n')
  f.write('\n')
  f.write('queue\n')

def print_condor_submit_file(f):
    f.seek(0)
    for line in f:
      print line,


if __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('evttype', help='evttype to cross validate. 1-5.', type=int)
  parser.add_argument('colnum', help='column index cross validate. 0 or 1.', type=int)
  parser.add_argument('--begin', type=float, default=0.1,
                               help='beginning bandwidth. ')
  parser.add_argument('--end', type=float, default=0.0001,
                               help='ending bandwidth. ')
  parser.add_argument('--step', type=float, default=0.0001,
                               help='decrement step size from start to end. ')
  parser.add_argument('--data_dir', default='kde_data', 
                                    help='Absolute path to kde input data.')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()

  dt = datetime.datetime.now()
  dirname = 'cv_output/evttype{0}_col{1}'.format(args.evttype, args.colnum)
  dirname = dirname + '/{0}{1}{2}_{3}{4}{5}'.format(dt.month, dt.day, dt.year, 
                                                    dt.hour, dt.minute, dt.second)
  os.makedirs(dirname) 

  f = tempfile.NamedTemporaryFile()
  write_condor_submit_file(f, dirname, args.data_dir, args.evttype, args.colnum, 
                           args.begin, args.end, args.step)
  if args.debug: 
    print_condor_submit_file(f)
  f.seek(0)
  subprocess.check_call(['condor_submit', f.name])
  f.close()
