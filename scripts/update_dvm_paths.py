#! /usr/bin/env python3

import os
import xml.etree.ElementTree as etree
import argparse
from distutils.util import strtobool

def main():
    parser = argparse.ArgumentParser('Update file paths in an DVM++ xml to the current directory')
    parser.add_argument('xml_file', help='the xml file to be updated', nargs='+')
    parser.add_argument('--output', help='alternative output directory')
    parser.add_argument('--input', help='alternative input directory')

    args = parser.parse_args()

    try:
        pwd = os.getcwd()

        for xml_file in args.xml_file:
            tree = etree.parse(xml_file)
            input_file = tree.find('io').find('input_dir')
            if args.input:
                input_file.set('string', args.input)
            else:
                input_file.set('string', pwd)

            output_file = tree.find('io').find('output_dir')
            if args.output:
                output_file.set('string', args.output)
            else:
                output_file.set('string', os.path.join(pwd, 'outputs'))

            tree.write(xml_file)

        outdir = output_file.get('string')
        if not os.path.isdir(outdir):
            sure = input('Output directory does not exist: ' + outdir + '. Create it? [y/n]   ')
            if strtobool(sure):
                os.mkdir('outputs')

        inputdir = tree.find('io').find('input_dir').get('string')
        domfile = tree.find('io').find('domain_file').get('string')
        dompath = os.path.join(inputdir, domfile)
        if not os.path.isfile(dompath):
            raise ValueError("Domain file not present: " + dompath)

    except (ValueError, etree.ParseError) as e:
        print("Error: " + str(e))

main()
