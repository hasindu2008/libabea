import sys
import os
import argparse
import numpy as np
import abea
import h5py


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
    NAME = "read2sigAlign"
    VERSION = "0.1.0"


    parser = MyParser(
        description="motif2sig - Sequence motif to signal")
    # group = parser.add_mutually_exclusive_group()
    parser.add_argument("-f", "--fastq",
                        help="input fastq file")
    parser.add_argument("-p", "--fast5",
                        help="fast5 file")
    parser.add_argument("-s", "--save",
                        help="save align file")
    parser.add_argument("--RNA", action="store_true",
                        help="Reads are from RNA sample")


    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)


    a = open(args.save, 'w')
    fastq_data = read_fastq(args.fastq)
    fast5_filepath = args.fast5
    f5_data = read_multi_fast5(fast5_filepath)
    for data in fastq_data:
        for readID in list(data.keys()):
            if readID in f5_data.keys():
                print(readID, data[readID]['seq'],
                      f5_data[readID]['digitisation'], f5_data[readID]['offset'],
                      f5_data[readID]['range'], f5_data[readID]['sampling_rate'], args.RNA, sep="\t")
                segs = get_segments(readID, data[readID]['seq'], f5_data[readID]['signal'],
                                    f5_data[readID]['digitisation'], f5_data[readID]['offset'],
                                    f5_data[readID]['range'], f5_data[readID]['sampling_rate'], args.RNA)
                # print(segs)
                # sys.exit()
                if segs is None:
                    # test for segs['FAIL'] too
                    print_err("Failed alignment: {}".format(readID))
                    continue
                else:
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
                segs = None


    a.close()

def get_segments(readID, seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE, RNA):

    ret = abea.abea_python(seq, SAMPLES, DIGITISATION, OFFSET, RANGE, SAMPLE_RATE, RNA=RNA)
    # print ("base_index start_raw_index(inclusive) end_raw_index(non_inclusive)")
    if 'FAIL' in ret:
        print_err("get_segments: failed to get alignment: {}".format(readID))
        # sys.exit()
        return None
    return ret


def read_fastq(filename, batch=1):
    '''
    read fastq
    '''
    c = 0
    b = 0
    keys = ['seq_length', 'seq']
    dic = {}
    with open(filename, 'r') as f:
        for ln in f:
            c += 1
            l = ln.strip('\n')
            if c == 1:
                b += 1
                s = l.split(' ')
                idx = s[0][1:]
                dic[idx] = {}
                dic[idx] = {k: None for k in keys}

            elif c == 2:
                    dic[idx]['seq'] = l
                    dic[idx]['seq_length'] = len(dic[idx]['seq'])

            if c >= 4:
                c = 0
                if b >= batch:
                    yield dic
                    b = 0
                    dic = {}



def read_multi_fast5(filename):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    hdf = h5py.File(filename, 'r')

    for read in list(hdf.keys()):
        readID = read.strip("read_")
        # read = "read_{}".format(readID)
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




if __name__ == '__main__':
    main()
