import sys
import os
import datetime
import tempfile
import subprocess

def write_condor_submit_file(
    f, data_dir, dirname, node_idx,
    first_idx, last_idx):

  result_fname = os.getcwd() + '/{0}/result/{1}.csv'.format(dirname, node_idx)
  out_fname = '{0}/out/{1}.out'.format(dirname, node_idx)
  err_fname = '{0}/err/{1}.err'.format(dirname, node_idx)
  log_fname = '{0}/log/{1}.log'.format(dirname, node_idx)
  data_dir = os.getcwd() + '/' + data_dir

  f.write('executable = cache_akde_score\n')
  f.write('universe = vanilla\n')
  f.write('arguments = {0} {1} {2} {3}\n'.format(data_dir, result_fname, first_idx, last_idx))
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
  parser.add_argument('--data_dir', default='kde_data', 
                                    help='Absolut path to kde input data.')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()

  dt = datetime.datetime.now()
  dirname = '{0}{1}{2}/{3}{4}{5}'.format(dt.month, dt.day, dt.year, 
                                         dt.hour, dt.minute, dt.second)
  for dir in [ 'result', 'out', 'err', 'log' ]:
    os.makedirs(dirname + '/' + dir) 

  total = 0
  with open('kde_data/test.dither.csv') as f:
    for line in f: total += 1

  per_node = total / args.nodes
  for i in range(args.nodes):
    output_fname = '{0}/result/{1}.csv'.format(dirname, i)
    first, last = i*per_node, (i+1)*per_node - 1
    if i == args.nodes - 1: last = total - 1

    f = tempfile.NamedTemporaryFile()
    write_condor_submit_file(f, args.data_dir, dirname, i, first, last)
    if args.debug and i == 0: 
      print_condor_submit_file(f)
    f.seek(0)

    subprocess.check_call(['condor_submit', f.name])
    f.close()
