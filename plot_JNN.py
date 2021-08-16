import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2020

    script description




    ----------------------------------------------------------------------------
    version 0.0 - initial



    TODO:
        - Take the whole array of segments, and overlay it on the raw signal
        - Makes for a good figure

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
    NAME = "plot_JNN"
    VERSION = "0.1.0"


    parser = MyParser(
        description="plot_JNN - for checking barcode segments")
    # group = parser.add_mutually_exclusive_group()
    parser.add_argument("-s", "--squigs",
                        help="squiggle file")
    parser.add_argument("-g", "--segs",
                        help="segment file")

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

    #main logic goes here

    seg_dic = get_segs(args.segs)
    c = 0
    with open(args.squigs, 'r') as f:
        for l in f:
            c += 1
            l = l.strip("\n")
            l = l.split("\t")
            readID = l[0]
            sig = np.array([float(i) for i in l[1:]], dtype=float)

            fig = plt.figure(1)
            #fig.subplots_adjust(hspace=0.1, wspace=0.01)
            ax = fig.add_subplot(111)
            length = len(sig)
            # Show segment lines
            for base, index, b_pose, r_start, r_end in seg_dic[readID]:
                i = int(r_start)
                j = int(r_end)
                if i > j:
                    i, j = j, i
                if i == j:
                    j = j + 1
                half_diff = (j-i) / 2
                mid = i + half_diff
                print(i, j, length)
                print(sig[i:j])
                y = max(sig[i:j]) + 20
                ax.text(mid, y, base)
                ax.axvline(x=i, color='m')
                ax.axvline(x=j, color='m')

            plt.plot(sig, color='k')
            plt.show()
            # plt.savefig("test_{}.png".format(c))
            plt.clf()
            sys.exit()



#def funct1():
    #blah blah
    #return

#...

def get_segs(filename):
    '''
    get segs and return dic with list items
    '''
    dic = {}
    with open(filename, 'r') as f:
        for l in f:
            l = l.strip('\n')
            l = l.split('\t')
            if l[0] not in list(dic.keys()):
                dic[l[0]] = []
            dic[l[0]].append(l[1:])
    return dic

if __name__ == '__main__':
    main()
