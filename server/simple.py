import paramiko

hostname = raw_input('Hostname: ')
user = raw_input('Username: ')
pwd = raw_input('Password: ')
client = paramiko.SSHClient()
client.load_system_host_keys()

client.connect(hostname, username=user, password=pwd)
stdin, stdout, stderr = client.exec_command('ls')
for line in stdout:
    print '... ' + line.strip('\n')

client.exec_command('python test.py')

client.close()
