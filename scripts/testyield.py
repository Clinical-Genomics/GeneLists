#!/usr/bin/env python
# encoding: utf-8

import sys

def printer():
    print('hey!')

def yielder():
    printer()
    for i in range(0,10):
        yield i

def main(argv):
    for yielded in yielder():
        print(yielded)

if __name__ == '__main__':
    main(sys.argv[1:])
