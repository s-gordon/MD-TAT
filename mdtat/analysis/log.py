#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import logging
import argparse
import sys


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def set_verbosity(parser, arg_verbosity):
    # Return argparse help if no arguments are specified. Catches situations
    # where only verbosity is set and does the same thing.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if arg_verbosity is True:
            parser.print_help()
            sys.exit(1)
    # Logger-control of feed-out
    # Normal output is prescribed to info
    # Debug output (accessed through -v/--verbose) is prescribed to debug
    """
    If verbosity set, change logging to debug.
    Else leave at info
    """
    if arg_verbosity is True:
        logging.basicConfig(format='%(levelname)s:\t%(message)s',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s:\t%(message)s',
                            level=logging.INFO)
