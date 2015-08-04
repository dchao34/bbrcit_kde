import sys
import os

if __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('result_dir', help='Directory containing the results.')
  parser.add_argument('output_fname', help='File to save collected results.')
  args = parser.parse_args()

  file_list = os.listdir(args.result_dir)
  with open(args.output_fname, 'w') as w:
    for fname in file_list:
      sys.stdout.write('Copying {0}...\n'.format(fname))
      sys.stdout.flush()
      with open(args.result_dir + '/' + fname, 'r') as f:
        for line in f:
          w.write(line)
  print 'Done.'
