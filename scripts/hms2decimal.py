#!/usr/bin/env python

import sys
import re

"""Handles conversion fo RA and Dec in various formats to decimal degrees"

# get file name
if len(sys.argv) == 1:
    ra  = raw_input('RA (hh:mm:ss.ss)    : ')
    dec = raw_input('Dec (+/-dd:mm:ss.ss): ')
elif len(sys.argv)== 2:
    ra = sys.argv[1]
    dec = raw_input('Dec (+/-dd:mm:ss.ss): ')
elif len(sys.argv)== 3:
    ra  = sys.argv[1]
    dec = sys.argv[2]
elif len(sys.argv)== 7:
    ra  = sys.argv[1] + ':' + sys.argv[2] + ':' + sys.argv[3]  
    dec = sys.argv[4] + ':' + sys.argv[5] + ':' + sys.argv[6]  
else:
    print 'usage: ra dec'
    exit(1)

# parse the RA

m = re.match('(\d\d)[:\s](\d\d)[:\s](\d\d(?:\.\d*)?)', ra)
if not m:
    print 'Could not interpret RA = ' +  ra
    print 'Example of correct format = 13:02:12.01'
    exit(1)

rad = 15.*(float(m.group(1)) + float(m.group(2))/60 + float(m.group(3))/3600)
print 'Decimal RA =  ' + str(rad)

# parse the Dec

m = re.match('([\+-])(\d\d)[:\s](\d\d)[:\s](\d\d(?:\.\d*)?)', dec)
if not m:
    print 'Could not interpret Dec = ' + dec
    print 'Example of correct format = -23:13:45.2'
    exit(1)

decd = float(m.group(2)) + float(m.group(3))/60 + float(m.group(4))/3600
if m.group(1) == '-':
    decd *= -1
    
print 'Decimal Dec = ' + str(decd)


