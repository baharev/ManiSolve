# Copyright (C) 2016, 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from datetime import datetime
from functools import partial, wraps
from glob import glob
from gzip import open as gz_open
from itertools import islice, tee
from os import makedirs, remove
from os.path import isdir
from pickle import load, dump, HIGHEST_PROTOCOL
from shutil import rmtree
import sys

from py3compat import izip, StringIO


def get_lines(filename):
    with open(filename, 'r') as f:
        return f.read().splitlines()

def as_pairs(itr):
    # [1, 2, 3, 4] -> (1, 2), (3, 4)
    return izip(islice(itr, 0, None, 2), islice(itr, 1, None, 2))

def pairwise(iterable):
    '''A generator object is returned.
    []  pairwise: []
    [1] pairwise: []
    [1,2,3] pairwise: [(1, 2), (2, 3)].'''
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def duplicates(iterable):
    seen = set()
    seen_add = seen.add
    return sorted(set(e for e in iterable if e in seen or seen_add(e)))

warning = partial(print, file=sys.stderr)

def print_timing(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        start = datetime.now() 
        print(start.strftime('%Y-%m-%d %H:%M:%S'))
        result = f(*args, **kwargs)
        end = datetime.now()
        print('Took: %.1f s' % (end-start).total_seconds())
        print(end.strftime('%Y-%m-%d %H:%M:%S'))
        return result
    return wrapper

#-------------------------------------------------------------------------------

class StdoutHijack:

    # If we do StdoutHijack again under a StdoutHijack with block in tee mode
    # then the output is repeated not once but twice, hence this hack: 
    __orig_stdout = sys.stdout

    def __init__(self):
        self.old = None
        self.log = StringIO()
    
    def __enter__(self):
        self.old = sys.stdout
        sys.stdout = self
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.old
    
    def write(self, message):
        if self.old is StdoutHijack.__orig_stdout:
            #self.old.write(message)  # uncomment if a tee like behavior is needed
            pass 
        self.log.write(message)
        
    def captured_text(self):
        return self.log.getvalue()

#-------------------------------------------------------------------------------

def serialize(obj, filename):
    with gz_open(filename, 'wb') as f:
        dump(obj, f, HIGHEST_PROTOCOL)

def deserialize(filename):
    with gz_open(filename, 'rb') as f:
        return load(f)

#-------------------------------------------------------------------------------

def clean(directory):
    if isdir(directory):
        print('Deleting folder "{}"'.format(directory))
        rmtree(directory)
    print('Creating folder "{}"'.format(directory))
    makedirs(directory)

def delete_files(pattern):
    for f in glob(pattern):
        remove(f)
