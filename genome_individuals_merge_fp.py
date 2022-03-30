import logging
import tarfile, io
import random
from storage.pyStorage import pyStorage
from Inspector import *

def lambda_handler(event, context):
    logging.info('Python HTTP trigger function processed a request.')
    inspector = Inspector()
    inspector.inspectAll()
    
    inspector.addAttribute("function_name", "genome_individual_merge_nv")
    inspector.addTimeStamp("Time Stamp at start")
 
    bucket_name = event["output_buckets"][0][3]
    
    key_columnsfile = event['key_columnsfile']
    
    pyStorage.create_cred_file(event['aws_access_key_id'], event['aws_secret_key'], event['aws_session_token'], event['gcp_client_email'], event['gcp_private_key'], event['gcp_project_id'])

 
    #Expected failure probability in percent
    failure_prob = event['failure_prob']

    #introduce failure probability
    if random.random() < (failure_prob/100):
        raise Exception('failure')


    key_datafiles = event['key_datafiles']
    
    c = 22

    columns = readColumns( key_columnsfile, inspector) 
    

    tars=[]
    for key in key_datafiles:
        tars.append(readTar(key, inspector))
    
    gnomes=[]
    for i in range(50):
        all_from_one_gnome = []
        for j in range(len(tars)):
            all_from_one_gnome.extend(tars[j][i])
        gnomes.append(all_from_one_gnome)
    
    print(len(gnomes), len(gnomes[0]))
    
    result_key = bucket_name + "tmp/chr"+str(c)+"n.tar.gz"
    save_concat_tar(result_key, gnomes, columns , c, inspector)
    
    inspector.addTimeStamp("Time Stamp at end")
    
    inspector.addAttribute("merged_result", result_key)
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
   
    
def readColumns(key_columnsfile, inspector):
    ##Read columns
    tmp_file = "/tmp/columns_file.txt"
    
    inspector.addTimeStamp("Time Stamp before Download 1")
    pyStorage.copy(key_columnsfile, tmp_file)
    inspector.addTimeStamp("Time Stamp after Download 1")
    
    object = open(tmp_file, "r").read()
    bytes(object, "utf-8")

    return object.rstrip().split("\t")

    
def readTar(key_gnomefile, inspector):
    print("---------", key_gnomefile)
    tmp_file = "/tmp/key_genome.tar.gz"
    
    inspector.addTimeStamp("Time Stamp before Download: " + key_gnomefile)
    pyStorage.copy(key_gnomefile, tmp_file)
    inspector.addTimeStamp("Time Stamp after Download: " + key_gnomefile)
    
    tar = tarfile.open(name= tmp_file, mode = "r:gz")
    files = []
    for member in tar.getmembers():
        f = tar.extractfile(member)
        if f == None:
            continue
        f= f.read().decode("utf-8").split("\n")
        files.append(f)
    if(not len(files)==50):
        print("Tar file has not enough files!")
    return files
    
def save_concat_tar(key, gnomes, columns, c, inspector):
    
    tmp_file = "/tmp/gnomes.tar.gz"
    
    file_out = io.BytesIO()
    tar = tarfile.open(name= tmp_file, mode = "w:gz")
    for i in range(len(gnomes)):
        gnomes[i]="\n".join(gnomes[i])
        string = io.BytesIO(gnomes[i].encode())
        info = tarfile.TarInfo(name="chr"+str(c)+"."+columns[i])
        info.size=len(string.getbuffer())
        tar.addfile(tarinfo=info, fileobj=string)
    tar.close()
    
    inspector.addTimeStamp("Time Stamp before Upload")
    pyStorage.copy(tmp_file, key)
    inspector.addTimeStamp("Time Stamp after Upload")
    
    return  key