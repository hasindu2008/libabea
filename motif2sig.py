import sys
import os
import argparse
import numpy as np
import abea

'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2020

    script description

    motif2sig.py takes the start and stop positions in a nucleotide sequence,
    and returns the raw signal co-ordinated associated with those nucleotides


    ----------------------------------------------------------------------------
    version 0.0 - initial



    TODO:
        -

    ----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2020 James Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''



# READ_LEN=7831 int
# N_SAMPLES=65030 int
# DIGITISATION=8192.0 float
# OFFSET=26.0 float
# RANGE=1444.9 float
# SAMPLE_RATE=4000.0 float

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def print_verbose(message, level):
    '''verbose printing'''
    if VERBOSE >= level:
        sys.stderr.write('info: %s\n' % message)

def print_err(message):
    '''error printing'''
    sys.stderr.write('error: %s\n' % message)

def main():
    '''
    do the thing
    '''
    NAME = "motif2sig"
    VERSION = "0.1.0"


    parser = MyParser(
        description="motif2sig - Sequence motif to signal")
    # group = parser.add_mutually_exclusive_group()
    parser.add_argument("-f", "--fastq",
                        help="input fastq file")
    parser.add_argument("-s", "--seq_sum",
                        help="input sequencing_summary file")
    parser.add_argument("-p", "--fast5",
                        help="fast5 file top path")

    parser.add_argument("-v", "--verbose", type=int, default=1,
                        help="Engage higher output verbosity")
    parser.add_argument("-V", "--version", action="store_true",
                        help="Engage higher output verbosity")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    global VERBOSE
    VERBOSE = 0
    if args.verbose:
        VERBOSE = args.verbose
        print_verbose("Verbose level {} active - dumping info to stderr".format(VERBOSE), 1)
        print_verbose("{} - {}".format(NAME, VERSION), 1)
        print_verbose("arg list: {}".format(args), 1)

    if args.version:
        sys.stderr.write("{} - {}\n".format(NAME, VERSION))
        sys.exit(1)


    # read seq_sum to map read -> fast5 (get pass state too for duplicates)
    reads2fast5 = read_seq_sum(args.seq_sum, args.fast5)
    # for readID, open fast5, and parse data to abea
    fastq_data = read_fastq(args.fastq)

    print("readID", "barcode", "cell", "score", "direction", "seq_length", "nt_start",
          "nt_stop", "raw_start", "raw_stop", sep="\t")
    for readID, data in list(fastq_data.keys):
        # cut out stop/start co-ords for raw signal, and dump into file

        fast5_filename = reads2fast5[readID][0]
        fast5_filepath = reads2fast5[readID][1]
        f5_data = read_multi_fast5(fast5_filepath)
        segs = get_segments(readID, data['seq'], f5_data['signal'],
                            f5_data['digitisation'], f5_data['offset'],
                            f5_data['range'], f5_data['sampling_rate'])

        nt_start = data['start']
        nt_stop = data['stop']
        raw_start = segs[nt_start][1]
        raw_stop = segs[nt_stop][2]

        seq_length = data['seq_length']
        barcode = data['barcode']
        cell = data['cell']
        score = data['score']
        direction = data['direction']

        # readID, nt_start, nt_stop, raw_start, raw_stop
        print(readID, barcode, cell, score, direction, seq_length, nt_start,
              nt_stop, raw_start, raw_stop, sep="\t")
        # plot seg cuts to visualise, I should be able to confirm with JNN



def read_seq_sum(filename, f5_path):
    '''
    read seq sum and pair read name with fast5 file name
    grab qscore field for pass/fail and legacy naming convention
    '''
    f5_dic = {}
    for dirpath, dirnames, files in os.walk(f5_path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)
                f5_dic[fast5] = fast5_file
    dic = {}
    qscore_col = False
    fast5_col = False
    readID_col = False
    head = True
    with open(filename, 'r') as f:
        for l in f:
            l = l.split('\t')
            if head:
                head = False
                for i in range(len(l)):
                    if l[i] == "mean_qscore_template":
                        qscore_col = i
                    if l[i] in ["filename_fast5", "filename"]:
                        fast5_col = i
                    if l[i] == "read_id":
                        readID_col = i

                if not qscore_col or not fast5_col or not readID_col:
                    print_err("columns not found: qscore: {} , fast5: {} , readID: {}". format(qscore_col, fast5_col, readID_col ))
                    sys.exit()
                continue
            dic[l[readID_col]] = [l[fast5_col], f5_dic[l[fast5_col]], float(l[qscore_col])]

    return dic


def read_fastq(filename, batch=1):
    '''
    read demuxed fastq file and yield generator with barcode, start, stop
    '''
    c = 0
    b = 0
    keys = ['barcode', 'cell', 'score', 'start', 'stop', 'direction', 'seq_length', 'seq']
    # dic = {'readID': None, 'barcode': None, 'start': None, 'stop': None,
    #        'score': None, 'seq_length': None, 'seq': 'None}
    dic = {}
    with open(filename, 'r') as f:
        for ln in f:
            c += 1
            l = ln.strip('\n')
            if c == 1:
                b += 1
                s = l.split(' ')
                idx = s[0][1:]
                k = [i.split('=') for i in s]
                dic[idx] = {keys: None for k in keys}
                for i in range(len(k)):
                    if k[i][0] in ["rev_bc", "fwd_bc"]:
                        dic[idx]['barcode'] = k[i][1]
                    elif k[i][0] in ["cellNum"]:
                        dic[idx]['cell'] = k[i][1]
                    elif k[i][0] in ["score"]:
                        dic[idx]['score'] = int(k[i][1])
                    elif k[i][0] in ["start"]:
                        dic[idx]['start'] = int(k[i][1])
                    elif k[i][0] in ["end"]:
                        dic[idx]['stop'] = int(k[i][1])
                    elif k[i][0] in ["rev_trimmed"]:
                        dic[idx]['seq'] = int(k[i][1])
                        dic[idx]['direction'] = "rev"
                    elif k[i][0] in ["fwd_trimmed"]:
                        dic[idx]['seq'] = int(k[i][1])
                        dic[idx]['direction'] = "fwd"


            elif c == 2:
                if dic[idx]['direction'] == "rev":
                    dic[idx]['seq'] = l + dic[idx]['seq']
                elif dic[idx]['direction'] == "fwd":
                    dic[idx]['seq'] = dic[idx]['seq'] + l
                else:
                    print_err("read_fastq: direction variable not valid: {}".format(direction))

                dic[idx]['seq_length'] = len(dic[idx]['seq'])
            if c >= 4:
                c = 0
                if b >= batch:
                    yield idx, dic
                    b = 0
                    dic = {keys: None for k in keys}



def read_multi_fast5(filename):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    with h5py.File(filename, 'r') as hdf:
        for read in list(hdf.keys()):
            f5_dic[read] = {'signal': [], 'readID': '', 'digitisation': 0.0,
                            'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}
            try:
                for col in hdf[read]['Raw/Signal'][()]:
                    f5_dic[read]['signal'].append(int(col))

                f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                f5_dic[read]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
                f5_dic[read]['offset'] = hdf[read]['channel_id'].attrs['offset']
                f5_dic[read]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
                f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                print_err("extract_fast5():failed to read readID: {}".format(read))
    return f5_dic

def get_segments(readID, seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE):

    ret = abea.abea_python(seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE)
    print ("base_index start_raw_index(inclusive) end_raw_index(non_inclusive)")
    if 'FAIL' in ret:
        print_err("get_segments: failed to get alignment: {}".format(readID))
        # sys.exit()
        return None

    return dic



if __name__ == '__main__':
    main()
