#!/software/containers/john_bercow.sif

# Code from Ryan to query IDT's website for sequence synthesiability.
import os
from pathlib import Path
import time
import json
import requests
from base64 import b64encode
from urllib import request, parse
import numpy as np

# ============================================
# FUNCTIONS
# ============================================
def vprint(str,verbose=False,**kwargs):
	if verbose: print(str,**kwargs)

def use_dir(dir):
	user_info_file = os.path.expanduser(os.path.join(dir, "info.json"))
	token_file = os.path.expanduser(os.path.join(dir, "token.json"))
	return user_info_file, token_file

def ask_for_user_data(user_info_file):
	user_info = {}
	idt_url = "https://www.idtdna.com/pages/tools/apidoc"
	msg = f"The first time you use the IDT " \
		f"complexity function requires you to " \
		f"input your IDT username, password, " \
		f"API client ID and API client secret. "\
		f"These will be stored securely in " \
		f"'{os.path.abspath(user_info_file)}', " \
		f"which only you have access to. Before " \
		f"you begin, go to {idt_url} and follow " \
		f"all the directions under the 'First time setup " \
		f"- Generating credentials' header. The client " \
		f"ID and client description can be " \
		f"anything. The client secret will be " \
		f"generated for you."

	print(msg)
	print("1) Please enter your IDT account username: ")
	user_info["username"] = input()
	print("2) Please enter your IDT account password: ")
	user_info["password"] = input()
	print("3) Please enter you API client ID (go to https://www.idtdna.com/pages/tools/apidoc): ")
	user_info["ID"] = input()
	print("4) Please enter your API client secret (go to https://www.idtdna.com/pages/tools/apidoc): ")
	user_info["secret"] = input()

	return user_info

def get_user_info(user_info_file):
	if os.path.exists(user_info_file):
		with open(user_info_file,'r') as f:
			user_info = json.load(f)

	else:
		#we have to set things up for the first time
		os.makedirs(os.path.dirname(user_info_file),exist_ok=True)

		user_info = ask_for_user_data(user_info_file)

		#do it this weird way just to make sure nobody can see your private stuff
		#simply creates an empty file (like touch)
		Path(user_info_file).touch()
		#make it so only the user can acces this file
		os.chmod(user_info_file,0o600)
		#now write the secret stuff ;)
		with open(user_info_file,'w+') as f:
			json.dump(user_info,f)

	return user_info

def delete_stored_token(token_file):
	if os.path.exists(token_file):
		os.remove(token_file)

def store_token(token,token_file):
	delete_stored_token(token_file)
	#do it this weird way just to make sure nobody can see your private stuff
	#creates an empty file
	Path(user_info_file).touch()
	#make it so only the user can acces this file
	os.chmod(user_info_file,0o600)
	#now write the secret stuff ;)
	with open(token_file,'w+') as f:
		json.dump(token,f)

# IDT's code
def get_new_token(user_info,verbose=False):
    """
    Create the HTTP request, transmit it, and then parse the response for the
    access token.

    The body_dict will also contain the fields "expires_in" that provides the
    time window the token is valid for (in seconds) and "token_type".
    """
    vprint("getting new token", verbose)
    client_id = user_info["ID"]
    client_secret = user_info["secret"]
    idt_username = user_info["username"]
    idt_password = user_info["password"]

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }

    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token",
                                    data = request_data,
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)
#    response = requests.get(post_request)
    # Process the HTTP response for the desired data
    body = response.read().decode()

    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)

    body_dict = json.loads(body)
    return body_dict

def get_stored_token(token_file):
	with open(token_file,'r') as f:
		return json.load(f)

def get_token(token_file, user_info,verbose=False):
	get_token_flag = False
	if os.path.exists(token_file):
                vprint(f"using token stored at {token_file}",verbose)
                modified_time = os.path.getmtime(token_file)
                current_time = time.time()
                token = get_stored_token(token_file)
                if current_time - modified_time > token["expires_in"]:
                        vprint("token expired",verbose)
                        get_token_flag = True
	else:
                vprint(f"no file found at {token_file}",verbose)
                get_token_flag = True

	if get_token_flag:
                token = get_new_token(user_info,verbose)
                vprint(f"storing token at {token_file}",verbose)
                store_token(token,token_file)
	return token

def query_complexity(seq, token,verbose=False):
	vprint("querying complexity",verbose)
	url = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenEBlockSequences"
	payload = f'[{{"Name":"My eBlock","Sequence":"{seq}"}}]'
	headers = {
		'Content-Type': 'application/json',
		'Authorization': f'Bearer {token}'
	}

	response = requests.request("POST", url, headers=headers, data=payload)

	return json.loads(response.text)

def total_score(dna_seq, idt_user_info):
    '''
    Tests synthesis complexity of a sequence.
    Returns the total score.
    '''
    idt_token = get_token(token_file, idt_user_info)
    response = query_complexity(dna_seq, idt_token["access_token"])
    
    if response == [[]]: # if sequence passes IDT's complexity test nothing is returned
        tot_score = 0.0

    else: # if not, returns a series of score (the sum total is what matters).
        tot_score = np.sum([x['Score'] for x in response[0]])

    return tot_score
