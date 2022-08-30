import os

for fname in os.listdir('.'):
    if 'shiftreg' not in fname:
        continue
    dst = 'memory_cell'.join(fname.split('shiftreg'))
    print('renaming: ', fname, ' -> ', dst)
    os.rename(fname, dst)