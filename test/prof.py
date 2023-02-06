#!/usr/bin/env python3

from pstats import Stats
import sys

s=Stats(sys.argv[1])
s.sort_stats('cumtime').print_stats(10)
