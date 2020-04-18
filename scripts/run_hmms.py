import subprocess
import sys, os, time, pysam
import argparse
import tempfile,shutil
from os.path import isfile, getsize

from ngs_utils import logger
from joblib import Parallel, delayed

def run_cmd(cmd):
    logger.info(cmd)
    subprocess.run(cmd, shell=True)

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--threads', type=int, action='store',
                      help='''Threads to use''', required=False, default = 1)
  parser.add_argument('-b', '--bamfile', type=str, action='store', 
                      help='''Bamfile to analyze''', required=True)                                              
  parser.add_argument('-d', '--directory', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  parser.add_argument('-H', '--hmm_list', default=None, type=str, required=True,
                      help='''File containing HMMs to use''')                      
  options = parser.parse_args()
  return options

def read_hmm_file(hmm_list):
  hmms = {}
  input = open(hmm_list, "r")
  idx = 0
  for line in input:
    #Todo, instead of index, note which viral family the hmm came from
    hmms[idx] = line.strip()
    idx+=1
  input.close()
  return hmms

def run_pipeline(options):
  hmms = read_hmm_file(options.hmm_list)
  if not os.path.exists(f'{options.directory}/logs'):
    os.mkdir(f'{options.directory}/logs')
  if not os.path.exists(f'{options.directory}/temp'):
    os.mkdir(f'{options.directory}/temp')
  if not os.path.exists(f'{options.directory}/temp/unmapped.fas'):
    print('Preparing unmapped sequences')
    prepare_unmapped_sequences(options)
  else:
    print(f'Unmapped sequences already prepared, skipping: {options.directory}/temp/unmapped.fas')
  #Now search HMMs against reads
  print("Running HMMs")
  start_time = time.time()

  threads_per_job = max(1, options.threads // len(hmms))
  def run_one_hmm(i):
    logger.info(f"\tRunning HMM {hmms[i]}" )
    output_file = f"{options.directory}/temp/hmmsearch.{i}"
    if not isfile(output_file) or getsize(output_file) == 0:
      run_cmd(f'nhmmer -o {output_file} --noali --cpu {threads_per_job} {hmms[i]} {options.directory}/temp/unmapped.fas')
    else:
      logger.info(f'\t{output_file} exists, skipping nhmer')
  Parallel(n_jobs=options.threads)(delayed(run_one_hmm)(i) for i in hmms.keys())

  end_time = time.time() - start_time
  print("Finished running against HMMs: %fs" % end_time)    
  print("Processing results\n")
  scores = {}
  for i in list(hmms.keys()):
    scores = read_nhmmer_result("%s/temp/hmmsearch.%s" % (options.directory, i), scores)
  output = open('%s/temp/reduced.csv' % options.directory, 'w')
  for score in list(scores.keys()):
    output.write("%s\n" % (','.join([str(s) for s in scores[score]])))
  output.close()

def read_nhmmer_result(file, scores):
  input = open(file, 'r')
  start_line = 'Query:'
  start = False
  for line in input:
    if start == False and line.find(start_line) == 0:
      start = True
      foo = next(input)
      foo = next(input)
      foo = next(input)      
    if start == True and line.find('read_') != -1:
      res = line.split()
      if scores.setdefault(res[3],[res[3], res[0], float(res[1])])[2] < float(res[1]):
        scores[res[3]] = [res[3], res[0], float(res[1])]
    elif line.strip() == "" and start == True:
      break
  return scores

def prepare_unmapped_sequences(options):
  start_time = time.time() 
  counter = 0
  bam = pysam.Samfile(options.bamfile, 'rb')
  fas = open('%s/temp/unmapped.fas' % options.directory, 'w')
  map = open('%s/temp/unmapped.map' % options.directory, 'w')
  
  for read in bam:
    if read.is_unmapped and not read.mate_is_unmapped and bam.references[read.rnext].find('chr') == 0:
      fas.write('>read_%d\n%s\n' % (counter, read.seq))
      map.write('%s\tread_%d\n' % (read.qname, counter))
      counter += 1
  fas.close()
  map.close()
  bam.close()
  end_time = time.time() - start_time
  print("Prepared sequences for searching against HMMs: %fs" % end_time)

if __name__ == '__main__': 
  start_time = time.time()
  reference_dir = os.environ['REFERENCE_REPO']
  options = parse_args()  
#   class Object(object):
#       pass
#   options = Object()  
#   options.bamfile='output.unknown.bam'
#   options.directory='tmp'
#   options.hmm_list='hmms.txt'
#   options.threads=1
  run_pipeline(options)