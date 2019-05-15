import os
import yaml
from subprocess import call
from utils import model_logger as logger


# user ssh to rivanna and submit job for user
PROJECT_DIR = os.path.dirname(os.path.realpath(__file__))

MONITOR_FILE = 'user_queue.yaml'
#USERCASE_DIR = '/home/yifan/Documents/Bart-web/BARTweb_docker/usercase' #path in the server
#LOG_FILE = '/home/yifan/Documents/Bart-web/BARTweb_docker/log/powhatan_trigger.log'

USERCASE_DIR = os.path.join(PROJECT_DIR, 'usercase')


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
                logger.info('Submit job: Job submission succeed!!!')
            else:
                logger.error('Submit job: Job submission failed!!!')
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
    commandret = 'I will write something when this actrually goes wrong'

    return callret,commandret



if __name__ == '__main__':
    # python trigger.py user_key
    import sys
    script_name = sys.argv[0]
    user_key = sys.argv[1]

    user_queue = load_queue() 
    check_queue(user_queue, user_key)    
    

