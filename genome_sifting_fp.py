import logging
import gzip
import random
from storage.pyStorage import pyStorage
from Inspector import *


def lambda_handler(event, context):

    logging.info('Python HTTP trigger function processed a request.')
    inspector = Inspector()
    inspector.inspectAll()

    inspector.addAttribute("function_name", "genome_sifting_nv")
    inspector.addTimeStamp("Time Stamp at start")
    
    bucket_name = event["output_buckets"][1]
    
    pyStorage.create_cred_file(event['aws_access_key_id'], event['aws_secret_key'], event['aws_session_token'], event['gcp_client_email'], event['gcp_private_key'], event['gcp_project_id'])


    #Expected failure probability in percent
    failure_prob = event['failure_prob']

    #introduce failure probability
    if random.random() < (failure_prob/100):
        raise Exception('failure')



    key_datafile = event['key'] #some tar
    
    c= 22

    content = get_data_file(key_datafile, inspector)
    
    
    print(content[0])

    nocomment = []
    lineNumber=0
    for x in content:
        if len(x)>0 and not x[0]=="#":
            lineNumber=lineNumber+1
            if "rs" in x and ("deleterious" in x or "tolerated" in x):
                nocomment.append([lineNumber] + x.rstrip().split("\t"))
                
    		
    
    for i in range(len(nocomment)):#len(nocomment)
    	info = nocomment[i][8].split("|")
    	try:	
    		sift = info[16].split("(")[1].split(")")[0]
    	except:
    		sift = ""
    	try:
    		phenotype = info[17].split("(")[1].split(")")[0]
    	except:
    		phenotype = ""
    	nocomment[i]=[nocomment[i][0],nocomment[i][3],info[4],sift,phenotype]

    print(nocomment[0])
    
    for i in range(len(nocomment)):
    	nocomment[i] = "\t".join(str(e) for e in nocomment[i])
    nocomment="\n".join(nocomment)
    
    content = nocomment.encode()
    
    content_file = open("/tmp/file.txt", "wb")
    content_file.write(content)
    content_file.close
    
    
    sifted = bucket_name+"/tmp/sifted.SIFT.chr"+str(c)+".txt"
    
    inspector.addTimeStamp("Time Stamp before Upload")
    pyStorage.copy("/tmp/file.txt", sifted)
    inspector.addTimeStamp("Time Stamp after Upload")
    
    
    inspector.addTimeStamp("Time Stamp at end")
    
    inspector.addAttribute("sifted", sifted)
    inspector.addAttribute("failure_prob", failure_prob)
    inspector.addAttribute("output_buckets", event['output_buckets'])
    inspector.addAttribute("aws_access_key_id", event['aws_access_key_id'])
    inspector.addAttribute("aws_secret_key", event['aws_secret_key'])
    inspector.addAttribute("aws_session_token", event['aws_session_token'])
    inspector.addAttribute("gcp_client_email", event['gcp_client_email'])
    inspector.addAttribute("gcp_private_key", event['gcp_private_key'])
    inspector.addAttribute("gcp_project_id", event['gcp_project_id'])
    inspector.inspectAllDeltas()
    
    return inspector.finish()
         


def get_data_file(key_datafile, inspector):
    tmp_file = "/tmp/file.vcf.gz"
    
    print("key_datafile: " + key_datafile)
    inspector.addTimeStamp("Time Stamp before Download")
    pyStorage.copy(key_datafile, tmp_file)
    inspector.addTimeStamp("Time Stamp after Download")
    
    with gzip.open(tmp_file, 'rb') as f:
        file_content = f.read().decode("utf-8") 

    return file_content.split('\n')
    