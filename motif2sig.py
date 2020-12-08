import sys
import os
import argparse
import numpy as np
import abea
import h5py
import traceback
# import matplotlib.pyplot as plt

'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2020

    script description

    motif2sig.py takes the start and stop positions in a nucleotide sequence,
    and returns the raw signal co-ordinated associated with those nucleotides


    ----------------------------------------------------------------------------
    version 0.0.0 - initial
    version 0.1.0 - scaling up for a full run



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
    parser.add_argument("-q", "--squiggle",
                        help="print raw squiggle integers to file")
    parser.add_argument("-g", "--sig_align",
                        help="print raw squiggle to base segments to file")
    parser.add_argument("-a", "--pA_convert", action="store_true",
                        help="convert squiggles to pA")
    parser.add_argument("-t", "--test", type=int, default=0,
                        help="number of reads to process when testing")

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
    if args.squiggle:
        w = open(args.squiggle, 'w')
    if args.sig_align:
        a = open(args.sig_align, 'w')
    count = 0
    for data in fastq_data:
        for readID in list(data.keys()):
            if args.test:
                if count >= args.test:
                    break
            # cut out stop/start co-ords for raw signal, and dump into file
            # print(data)
            fast5_filename = reads2fast5[readID][0]
            fast5_filepath = reads2fast5[readID][1]
            f5_data = read_multi_fast5(fast5_filepath, readID)
            segs = get_segments(readID, data[readID]['seq'], f5_data[readID]['signal'],
                                f5_data[readID]['digitisation'], f5_data[readID]['offset'],
                                f5_data[readID]['range'], f5_data[readID]['sampling_rate'])
            if segs is None:
                # test for segs['FAIL'] too
                print_err("Failed alignment: {}".format(readID))
                continue
            try:
                # index is off by 2 starting at 0
                nt_start = data[readID]['start'] - 2
                nt_stop = data[readID]['stop'] - 2
                if segs[nt_start][1] < 0:
                    raw_start = segs[nt_start-1][2]
                else:
                    raw_start = segs[nt_start][1]
                if segs[nt_stop][2] < 0:
                    raw_stop = segs[nt_stop+1][1]
                else:
                    raw_stop = segs[nt_stop][2]

                seq_length = data[readID]['seq_length']
                barcode = data[readID]['barcode']
                cell = data[readID]['cell']
                score = data[readID]['score']
                direction = data[readID]['direction']
            except:
                traceback.print_exc()
                print_err("Failed assignment: {}".format(readID))
                continue

            # readID, nt_start, nt_stop, raw_start, raw_stop
            print(readID, barcode, cell, score, direction, seq_length, data[readID]['start'],
                  data[readID]['stop'], raw_start, raw_stop, sep="\t")

            if args.squiggle:
                ar = []
                if args.pA_convert:
                    pa_sig = convert_to_pA(f5_data[readID])
                    for i in pa_sig:
                        ar.append(str(i))
                else:
                    for i in f5_data[readID]['signal']:
                        ar.append(str(i))
                w.write("{}\t{}\n".format(readID, '\t'.join(ar)))

            if args.sig_align:
                a.write("{}\t{}\t{}\t{}\t{}\t{}\n".format("readID", "base", "seg_index", "base_position", "sig_start", "sig_stop"))
                # get the positions then sort them to loop over
                segs_pos = list(segs.keys())
                segs_pos.sort()
                for pos in segs_pos:
                    i, j, k = segs[pos]
                    base = data[readID]['seq'][i]
                    try:
                        if segs[pos][1] < 0:
                            sig_start = segs[pos-1][2]
                        else:
                            sig_start = segs[pos][1]
                        if segs[pos][2] < 0:
                            sig_stop = segs[pos+1][1]
                        else:
                            sig_stop = segs[pos][2]
                        a.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(readID, base, pos, i, sig_start, sig_stop))
                    except:
                        print_err("{} pos missing: {}".format(readID, pos))
                        continue


            # plot seg cuts to visualise, I should be able to confirm with JNN
            fast5_filename = None
            fast5_filepath = None
            f5_data = None
            segs = None
            nt_start = None
            nt_stop = None
            raw_start = None
            raw_stop = None
            seq_length = None
            barcode = None
            cell = None
            score = None
            direction = None
        if args.test:
            if count >= args.test:
                break

        count += 1
    if args.squiggle:
        w.close()
    if args.sig_align:
        a.close()


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
    qscore_col = None
    fast5_col = None
    readID_col = None
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

                if qscore_col is None or fast5_col is None or readID_col is None:
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
            # print_err("read_fastq: c={}".format(c))
            # print_err("read_fastq: l={}".format(l))
            if c == 1:
                b += 1
                s = l.split(' ')
                idx = s[0][1:]
                k = [i.split('=') for i in s]
                dic[idx] = {}
                dic[idx] = {k: None for k in keys}
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
                        dic[idx]['seq'] = k[i][1]
                        dic[idx]['direction'] = "rev"
                    elif k[i][0] in ["fwd_trimmed"]:
                        dic[idx]['seq'] = k[i][1]
                        dic[idx]['direction'] = "fwd"


            elif c == 2:
                if dic[idx]['direction'] == "rev":
                    dic[idx]['seq'] = l + dic[idx]['seq']
                elif dic[idx]['direction'] == "fwd":
                    dic[idx]['seq'] = dic[idx]['seq'] + l
                else:
                    print_err("read_fastq: direction variable not valid: {} - {}".format(idx, dic[idx]['direction']))
                    # print_err("read_fastq: {}".format(dic))
                    # print_err("read_fastq: c={}".format(c))
                    # print_err("read_fastq: b={}".format(b))
                    # print_err("read_fastq: idx={}".format(idx))
                    # print_err("read_fastq: ln={}".format(ln))


                dic[idx]['seq_length'] = len(dic[idx]['seq'])

            if c >= 4:
                c = 0
                if b >= batch:
                    yield dic
                    b = 0
                    dic = {}



def read_multi_fast5(filename, readID):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    hdf = h5py.File(filename, 'r')

    # for read in list(hdf.keys()):
    read = "read_{}".format(readID)
    f5_dic[readID] = {'signal': [], 'readID': '', 'digitisation': 0.0,
                    'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}
    try:
        f5_dic[readID]['signal'] = hdf[read]['Raw/Signal'][()]
    except KeyError:
        hdf.close()
        filename = filename.replace("pass", "fail")
        hdf = h5py.File(filename, 'r')
    try:
        # for col in hdf[read]['Raw/Signal'][()]:
            # f5_dic[readID]['signal'].append(int(col))
        f5_dic[readID]['signal'] = hdf[read]['Raw/Signal'][()]

        f5_dic[readID]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
        f5_dic[readID]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
        f5_dic[readID]['offset'] = hdf[read]['channel_id'].attrs['offset']
        f5_dic[readID]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
        f5_dic[readID]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']

    except:
        hdf.close()
        traceback.print_exc()
        print_err("extract_fast5():failed to read readID: {}".format(read))

    hdf.close()

    return f5_dic

def get_segments(readID, seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE):

    ret = abea.abea_python(seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE)
    # print ("base_index start_raw_index(inclusive) end_raw_index(non_inclusive)")
    if 'FAIL' in ret:
        print_err("get_segments: failed to get alignment: {}".format(readID))
        # sys.exit()
        return None

    return ret

def convert_to_pA(d):
    '''
    convert raw signal data to pA using digitisation, offset, and range
    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        rawptr[j] = (rawptr[j] + offset) * raw_unit;
    }
    '''
    digitisation = d['digitisation']
    range = d['range']
    offset = d['offset']
    raw_unit = range / digitisation
    new_raw = []
    for i in d['signal']:
        j = (i + offset) * raw_unit
        new_raw.append("{0:.2f}".format(round(j,2)))
    return new_raw


if __name__ == '__main__':
    main()
