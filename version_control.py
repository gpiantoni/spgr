from datetime import datetime
from socket import gethostname
from subprocess import check_output
import phypno

phypno_path = phypno.__path__[0] + '/../.git'
phypno_ver = check_output('git --git-dir ' + phypno_path + ' rev-parse HEAD',
                          shell=True).decode('utf-8')

print('Last run on ' + gethostname() + ' at ' + str(datetime.now()) + '\n')
print('Phypno Version: ' + phypno_ver)
print(check_output('pip freeze', shell=True).decode('utf-8').replace('\n', ', '))
