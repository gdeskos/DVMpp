#! /usr/bin/env python3

import os
import xml.etree.ElementTree as etree
import argparse

def getElem(tree, elem):
    root = tree.getroot();
    for e in root.iter():
        if e.tag == elem:
            return e

    # only get here if we haven't found anything
    raise ValueError('cannot find field ' + elem)

def updateElem(elem, value):
    key = elem.keys()[0] # we know there is only ever one
    elem.set(key, value)


def main():
    parser = argparse.ArgumentParser('Update a field in an ICRD3D xml file')
    parser.add_argument('xml_file', help='the xml file to be updated', nargs='+')
    parser.add_argument('field', help='field to be updated')
    parser.add_argument('value', help='new value')
    parser.add_argument('-o', '--output', help='optional output file')

    args = parser.parse_args()

    try:
        for f in args.xml_file:
            tree = etree.parse(f)
            elem = getElem(tree, args.field)
            updateElem(elem, args.value)

            if args.output:
                tree.write(args.output)
            else:
                tree.write(f)

    except (IOError, ValueError) as e:
        print('Error: ' + str(e))

main()


