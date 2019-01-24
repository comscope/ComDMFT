import os, subprocess


def get_env_dict(key="SLURM_"):
    '''get a list of environment variables containing key word of key.
    '''
    envs = {}
    for e in os.environ:
        if key in e:
            envs[e] = os.environ[e]
    return envs


def unset_environ(envs):
    '''unset the list of environment variables envs.
    '''
    for e in envs:
        if e in os.environ:
            del os.environ[e]


def set_environ(envs):
    '''set the list of environment variables envs.
    '''
    for e in envs:
        if e not in os.environ:
            os.environ[e] = envs[e]

