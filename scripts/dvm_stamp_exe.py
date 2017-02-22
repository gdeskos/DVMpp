#! /usr/bin/env python3
'''Get the timestamp and execution time from the dvm screen output'''

import argparse


def main():
    parser = argparse.ArgumentParser(
        description='Get timestamp and execution time from the dvm screen'
        'output (eg DVM++ file.xml > file_out.txt)')
    parser.add_argument('files', help='files to process', nargs='+')
    parser.add_argument(
        '--matlab',
        help='print time stamps as cell array for matlab',
        action='store_true')

    args = parser.parse_args()

    xmls = []
    stamps = []
    rts = []
    for file in args.files:
        # Read the first 5 lines to get the xml filename and the stamp
        with open(file) as f:
            lines = [f.readline() for l in range(5)]

        xmls.append(lines[3].split()[3])
        stamps.append(lines[4].split()[3])

        # The run time is printed almost at the end of the file. If we have lots of
        # time steps the file could be massive so just seek straight to end and
        # come back a bit. We can only seek in binary mode so have to play some
        # games getting a regular string version back.
        with open(file, 'rb') as fb:
            fb.seek(0, 2)
            # 100 bytes arbitrary here, it is enough to capture the runtime line
            fb.seek(-100, 2)
            rts.append(
                fb.read(100).split(b'\n')[3].decode('ascii').split()[2]
                .replace('s', ''))

    longest = max(len(fname) for fname in args.files)
    print('File'.ljust(longest), 'Timestamp'.ljust(20), 'Execution time')
    for x, s, r in zip(xmls, stamps, rts):
        print(x.ljust(longest + 1), s.ljust(21), r, sep='')

    if args.matlab:
        print('')
        times = ["\'" + t + "\'" for t in stamps]
        print('{', ', '.join(times), '};', sep='')

        print('[', ', '.join(rts), ']', sep='')


main()
