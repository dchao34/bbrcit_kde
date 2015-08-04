import sys
import os
import datetime
import tempfile
import subprocess

def write_condor_submit_file(f, dirname, node_idx, cached_kde_fname, bags):

  result_fname = os.getcwd() + '/{0}/result/{1}.csv'.format(dirname, node_idx)
  out_fname = '{0}/out/{1}.out'.format(dirname, node_idx)
  err_fname = '{0}/err/{1}.err'.format(dirname, node_idx)
  log_fname = '{0}/log/{1}.log'.format(dirname, node_idx)
  cached_kde_fname = os.getcwd() + '/' + cached_kde_fname

  f.write('executable = mle_cvxnl.py\n')
  f.write('universe = vanilla\n')
  f.write('arguments = {0} {1} {2}\n'.format(bags, cached_kde_fname, result_fname))
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
  parser.add_argument('nodes', help='Number of nodes to use.', type=int)
  parser.add_argument('bags', help='Number of bag samples per node.', type=int)
  parser.add_argument('cached_kde_fname', 
                              help='Absolute path to cached kde scores.')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()

  dt = datetime.datetime.now()
  dirname = 'mle_output/{0}{1}{2}/{3}{4}{5}'.format(dt.month, dt.day, dt.year, 
                                                   dt.hour, dt.minute, dt.second)
  for dir in [ 'result', 'out', 'err', 'log' ]:
    os.makedirs(dirname + '/' + dir) 

  for i in range(args.nodes):
    f = tempfile.NamedTemporaryFile()
    write_condor_submit_file(f, dirname, i, args.cached_kde_fname, args.bags)
    if args.debug and i == 0: 
      print_condor_submit_file(f)
    f.seek(0)

    subprocess.check_call(['condor_submit', f.name])
    f.close()
