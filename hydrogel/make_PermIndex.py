import sys

args = sys.argv
index = int(args[1])

Np = 351000
Nt = 1821000

outf = 'PermIndex'
with open(outf, 'wt') as of:
    for i in range(1, Np+1):
        of.write('{0:8d}\t{1:8d}\n'.format(i,i))

    for i in range(Np+1, Np+4):
        of.write('{0:8d}\t{1:8d}\n'.format(i,i+3*(index-1)))

    for i in range(Np+4, Np+4+3*(index-2)):
        of.write('{0:8d}\t{1:8d}\n'.format(i,i))

    of.write('{0:8d}\t{1:8d}\n'.format(Np+4+3*(index-2),  Np+1))
    of.write('{0:8d}\t{1:8d}\n'.format(Np+4+3*(index-2)+1,Np+2))
    of.write('{0:8d}\t{1:8d}\n'.format(Np+4+3*(index-2)+2,Np+3))

    for i in range(Np+4+3*(index-2)+3,Nt+1):
        of.write('{0:8d}\t{1:8d}\n'.format(i,i))
