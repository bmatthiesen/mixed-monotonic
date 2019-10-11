# Copyright (C) 2018-2019 Bho Matthiesen, Christoph Hellings
# 
# This program is used in the article:
#
# Bho Matthiesen, Christoph Hellings, Eduard A. Jorswieck, and Wolfgang
# Utschick, "Mixed Monotonic Programming for Fast Global Optimization,"
# submitted to IEEE  Transactions on Signal Processing.
# 
# 
# License:
# This program is licensed under the GPLv2 license. If you in any way use this
# code for research that results in publications, please cite our original
# article listed above.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


import os
import os.path
import time
import random

def acquire(fn):
    while True:
        if not os.path.exists(fn):
            f = open(fn,'x')
            f.write(str(os.getpid()))
            f.close()
            break

        time.sleep(random.randint(1,100)/1000)

def release(fn):
    if os.path.exists(fn):
        f = open(fn,'r')
        if f.read() == str(os.getpid()):
            os.unlink(fn)
            return True
        else:
            return False
    else:
        return True

class Lock:
    def __init__(self, fn):
        self.__fn = fn

    def __enter__(self):
        acquire(self.__fn)

    def __exit__(self, *args, **kwargs):
        release(self.__fn)
