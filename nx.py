# Copyright (C) 2018 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>

# import sys
# from socket import gethostname
# sys.path.insert(0,  '/home/ali/nx/networkx' if gethostname() == 'X230' else '/users/ali/networkx')
# del sys 
# del gethostname

import sys
import os
sys.path.insert(0,  os.path.abspath('nx/networkx'))
del sys 
del os

from networkx import *
__version__ = networkx.__version__
