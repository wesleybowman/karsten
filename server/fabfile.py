from fabric.api import run, env, settings
from fabric.tasks import execute

def taskA():
    env.password = 'jgalt4u'
    #env.always_use_pty=False
    t = run('ls')
    return t

def taskB():
    env.password = 'jgalt4u'

    a = run('python fabTest.py', pty=True)
    print a
    return a

if __name__ == '__main__':
    t = execute(taskA, hosts=['107002b@julian.acadiau.ca'])
    a = execute(taskB, hosts=['107002b@julian.acadiau.ca'])
    print t
    print a
    print a['107002b@julian.acadiau.ca'].shape
    #taskA()
    #a = taskB()
