import os
import yaml
import logging
from subprocess import call

# user ssh to rivanna and submit job for user

MONITOR_FILE = 'user_queue.yaml'
USERCASE_DIR = '/home/yifan/Documents/Bart-web/BARTweb_docker/usercase' #path in the server
LOG_FILE = '/home/yifan/Documents/Bart-web/BARTweb_docker/log/powhatan_trigger.log'

# config logging
# classical logger
logger = logging.getLogger('trigger')
logger.setLevel(logging.DEBUG)
# write log to file
file_handler = logging.FileHandler(LOG_FILE)
file_handler.setLevel(logging.DEBUG)
# write log to console
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.ERROR)
# log formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

def load_queue():
    try:
        # load user queue file
        with open(os.path.join(USERCASE_DIR, MONITOR_FILE), 'r') as fopen:
            data = yaml.load(fopen)
    # handle exception
    except yaml.YAMLError as exc:
        logging.error("Error while parsing YAML file:")
        if hasattr(exc, 'problem_mark'):
            if exc.context != None:
                logger.error('  parser says\n' + str(exc.problem_mark) + '\n  ' +
                       str(exc.problem) + ' ' + str(exc.context) +
                       '\nPlease check data and retry.')
            else:
                logger.error('  parser says\n' + str(exc.problem_mark) + '\n  ' +
                       str(exc.problem) + '\nPlease correct data and retry.')
        else:
                logger.error("Something went wrong while parsing yaml file")
        return
    return data

def check_queue(user_queue, user_key):
    # queue containing the users' jobs which are left undone
    if user_key not in user_queue:
        logger.error('Fatal error: user_queue.yaml does not have {} in the list..'.format(user_key))
        return
    user_path = os.path.join(USERCASE_DIR, user_key)
    if os.path.exists(user_path):
        logger.info('Check user queue: ' + str(user_key) + ' directory exists.')
        logger.info('Submit job: Begin execute task for this user.')
        # start marge-bart pipeline according to config
        # if cleaned upload files
        user_config_file = os.path.join(user_path, 'user.config')
        user_config_data = {}
        with open(user_config_file, 'r') as fopen:
            user_config_data = yaml.load(fopen)
        if 'status' in user_config_data:
            return

        if user_queue[user_key]['status'] == 'Cleaned':
            callret,commandret = execute_sh(user_config_data)
            if callret == 0:
                logger.info('Submit job: SSH Rivanna succeed!!!')
            else:
                logger.error('Submit job: SSH Rivanna failed!!!')
                logger.error('error message from subprocess.call: {}'.format(callret))
                logger.error('error message from command line: {}'.format(commandret))
        else:
            logger.error('This user uploaded harmful stuff: {}'.format(user_key))
    else:
        logger.error('Check user queue: ' + str(user_key) + ' user file is missing. Please check!')

def execute_sh(user_data):

    # TODO: get the return message from the command line
    command_dir = os.path.join(user_data['user_path'],'run_bart.sh')
    cmd = '/bin/bash '+command_dir
    callret = call(cmd, shell=True)
    commandret = []

    return callret,commandret




def ssh_rivanna(user_key, user_path):
    # the slurm is generated when user finished configuration, and move with usercase
    user_config_file = os.path.join(user_path, 'user.config')

    rivanna_usercase_path = ''
    user_config_data = {}
    with open(user_config_file, 'r') as fopen:
        user_config_data = yaml.load(fopen)
    if 'status' in user_config_data:
        return 1 # if already processing, skip it

    rivanna_usercase_path = user_config_data['user_path']  # get path in Rivanna

    if not rivanna_usercase_path or not user_config_data:
        logger.error('SSH Rivanna: Rivanna user path does not exist, please check!')
        logger.error('SSH Rivanna: ' + user_key)
        if user_config_data:
            user_config_data['status'] = 'Error'

            with open(user_config_file, 'w') as fopen:
                yaml.safe_dump(user_config_data, fopen, default_flow_style=False, encoding='utf-8', allow_unicode=True)
        return -1    
    try: 
        # connect to rivanna
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.load_system_host_keys()
        client.connect(HOST, username=USER)
        logger.info('SSH Rivanna: Connecting to rivanna.')
        cmd = 'cd {0} && sbatch {0}/exe.slurm'.format(rivanna_usercase_path)
        logger.info('SSH Rivanna: ' + cmd)
        stdin, stdout, stderr = client.exec_command(cmd)
        logger.info(stdout.read())
        logger.info(stderr.read())

        client.close()

        user_config_data['status'] = 'Processing'
        with open(user_config_file, 'w') as fopen:
            yaml.safe_dump(user_config_data, fopen, default_flow_style=False, encoding='utf-8', allow_unicode=True)
        return 0

    except paramiko.BadHostKeyException as e:
        logger.error('SSH failed!!! BadHostKeyException')
        logger.error(e)
        raise Exception('BadHostKeyException')
    except paramiko.AuthenticationException as e:
        logger.error('SSH failed!!! AuthenticationException')
        logger.error(e)
        raise Exception('AuthenticationException')
    except paramiko.SSHException as e:
        logger.error('SSH failed!!! SSHException')
        logger.error(e)
        raise Exception('SSHException')
    except socket.error as e:
        logger.error('SSH failed!!! Communication socket error')
        logger.error(e)

    # (paramiko.BadHostKeyException, paramiko.AuthenticationException, 
    #                 paramiko.SSHException, socket.error) as e:
    #     logger.error('SSH failed!!! please check authorizations or internet connection:')
    #     logger.error(e)


if __name__ == '__main__':
    # python trigger.py user_key
    import sys
    script_name = sys.argv[0]
    user_key = sys.argv[1]

    user_queue = load_queue() 
    check_queue(user_queue, user_key)    
    

