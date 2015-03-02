#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import, unicode_literals
from datetime import datetime
import subprocess
import os

def getgittag(filename, date=None):
    """Gets the current version of a gene list. If date is provided,
    get commit on or before that date.

    Args:
        filename (str): the name of the gene list
        date (date, optional): a git parsable date, e.g. Mon Feb 9 14:19:16 2015 or 2015-02-09

    Returns (str): a version (tag) of the gene list

    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(filename))

    try:
        command = ['git', 'describe']
        if date:
            commit_sha1 = subprocess.check_output(['git', 'rev-list', '-n', '1', '--before=%s' % date, 'master']).decode('utf-8').strip()
            command.extend( ('--tags', commit_sha1) )
        tag = subprocess.check_output(command).decode('utf-8').strip()
    except subprocess.CalledProcessError:
        return None

    os.chdir(cwd)

    return tag

def getgitlastmoddate(filename, date_format='%Y%m%d'):
    """Gets the last modifiation date of a gene list

    Args:
        filename (str): the name of the gene list

    Returns (str): return date (e.g. 20150225)

    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(filename))
    full_str_date = subprocess.check_output(['git', 'log', '-1', '--format=%ad', '--', filename]).decode('utf-8').strip()
    os.chdir(cwd)

    # Mon Feb 9 14:19:16 2015 +0100
    full_date = datetime.strptime(full_str_date.partition('+')[0], '%a %b %d %H:%M:%S %Y ')

    return full_date.strftime(date_format)

def main(args):
    print(getgitlastmoddate(__file__, '%c'))

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
