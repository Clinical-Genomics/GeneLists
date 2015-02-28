#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import, unicode_literals
from datetime import datetime
import subprocess
import os

def getgittag(filename):
    """Gets the current version of a gene list

    Args:
        filename (str): the name of the gene list

    Returns (str): a version (tag) of the gene list

    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(filename))

    try:
        tag = subprocess.check_output(['git', 'describe']).decode('utf-8').strip()
    except subprocess.CalledProcessError:
        return None

    os.chdir(cwd)

    return tag

def getgitlastmoddate(filename):
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

    return full_date.strftime('%Y%m%d')
