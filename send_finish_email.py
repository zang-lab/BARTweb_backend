# for sending e-mail to user
import os, sys
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import yaml

def get_config(user_path):
    '''
    load data from user.config

    '''
    config_file = os.path.join(user_path, 'user.config')

    # handle non-existing exception
    if not os.path.exists(config_file):
        return None

    with open(config_file, 'r') as fopen:
        user_data = yaml.load(fopen)

    return user_data

def send_email(user_email,title,message):
	HOST_ADDRESS = "zanglab.service@gmail.com"
	PASSWORD = "ZangLab2018"
	msg = MIMEMultipart()
	msg['From'] = HOST_ADDRESS
	msg['To'] = user_email
	msg['Subject'] = title
	msg.attach(MIMEText(message, 'plain'))
	server = smtplib.SMTP_SSL("smtp.gmail.com")
	try:
	    server.login(HOST_ADDRESS, PASSWORD)
	    msg = msg.as_string()
	    server.sendmail(HOST_ADDRESS, user_email, msg)
	except smtplib.SMTPAuthenticationError:
	    return False, "username or password is wrong"
	except:
		return False, "errors in sending key to e-mail..."
	finally:
	    server.quit()  # finally close the connection with server

	return True, "send e-mail to user successfully..."



def if_success(user_path,dataType):
	done = False
	bart_output_dir = user_path+'/download/'
	files = os.listdir(bart_output_dir)
	count = 0
	if dataType == 'Geneset':
	    for file in files:
	        if '_auc.txt' in file: 
	            count = count+1
	        if '_bart_results.txt' in file: 
	            count = count+1
	        if '_adaptive_lasso_Info.txt' in file: 
	            count = count+1
	        if '_enhancer_prediction_lasso.txt' in file: 
	            count = count+1
	    if count==4:
	        done = True
	if dataType == 'ChIP-seq':
	    for file in files:
	        if '_auc.txt' in file: 
	            count = count+1
	        if '_bart_results.txt' in file: 
	            count = count+1
	    if count==2:
	        done = True
	return done


if __name__ == '__main__':
	path = sys.argv[1]
	user_data = get_config(path)
	key = user_data['user_key']
	email = user_data['user_email']
	dataType = user_data['dataType']
	successful_message = '''
success {}
'''.format(key)
	unsuccessful_message = '''
no success {}
'''.format(key)
	if email != '':
		if if_success(path,dataType):
			title = 'BART done'
			result,message = send_email(email,title,successful_message)
			print(message)
		else:
			title = 'BART error'
			result,message = send_email(email,title,unsuccessful_message)
			print(message)