#!/usr/bin/python

import os
import sys
import random

def main(n, p):

    for i in range(int(p)):
        with open("/tmp/fredericmarchand/input-" + str(i) + ".txt", 'w') as f:

            f.write(n + '\n')
            f.write(p + '\n')

            for j in range(int(n)/int(p)):
                r = random.randint(1, 100000)
                f.write(str(r) + '\n')
        
            f.close()

if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print "Wrong number of arguments"
