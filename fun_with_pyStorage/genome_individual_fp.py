import logging
import tarfile
import io
import gzip
import random
from Inspector import *
from storage.pyStorage import pyStorage

def lambda_handler(event, context):
    
    logging.info('Python HTTP trigger function processed a request.')
    
    inspector = Inspector()
    inspector.inspectAll()
    
    inspector.addTimeStamp("Time Stamp at start")
    
    pyStorage.create_cred_file(event['aws_access_key_id'], event['aws_secret_key'], event['aws_session_token'], event['gcp_client_email'], event['gcp_private_key'], event['gcp_project_id'])
    
    
    key_datafile = event['key_input']

    k = event['k']

    x = event['x']
    
    #Expected failure probability in percent
    failure_prob = event['failure_prob']
    
    #introduce failure probability
    if random.random() < (failure_prob/100):
        raise Exception('failure')

    
    bucket_name = event["output_buckets"][2]
    
    key_columnsfile = event['key_columnsfile']
    
    
    counter = int((x*10000/k)+1)
    
    stop = int((x+1)*(10000/k)+1)

    c = 22


    #try:
    columns = readColumns(key_columnsfile, inspector)
    
    lines = readGzip(key_datafile, inspector)
    
    lines = lines[counter-1:stop-1]
    nocomment = []
    for x in lines:
    	if not x[0]=="#":
    		nocomment.append(x.rstrip().split("\t"))

    gnomes = doIndividual(nocomment)

    key_result = bucket_name+"tmp/chr"+str(c)+"n-"+str(counter)+"-"+str(stop)+".tar.gz"    #chr21n-1-1001.tar.gz
    
    tarfile_name = createTarWithGnomes(gnomes,columns,c)
    
    
    inspector.addTimeStamp("Time Stamp before Upload")
    
    pyStorage.copy(tarfile_name, key_result)
    
    inspector.addTimeStamp("Time Stamp after Upload")
    
    inspector.addAttribute("key_columnsfile", key_columnsfile)
    inspector.addAttribute("key_result", key_result)
    inspector.addAttribute("failure_prob", failure_prob)
    inspector.addAttribute("output_buckets", event['output_buckets'])
    inspector.addAttribute("aws_access_key_id", event['aws_access_key_id'])
    inspector.addAttribute("aws_secret_key", event['aws_secret_key'])
    inspector.addAttribute("aws_session_token", event['aws_session_token'])
    inspector.addAttribute("gcp_client_email", event['gcp_client_email'])
    inspector.addAttribute("gcp_private_key", event['gcp_private_key'])
    inspector.addAttribute("gcp_project_id", event['gcp_project_id'])
    
    
    inspector.addTimeStamp("Time Stamp at end")
    
    
    inspector.inspectAllDeltas()
    

    return inspector.finish() 

def doIndividual(nocomment):
    individuals=0
    gnomes=[]
    for j in range(9,59):#2513
        gnome=[]
        for i in range(len(nocomment)):
            e = nocomment[i]
            cs = e[j].split("|")
            val = e[7].split(";")[8].split("=")[1]
            	
            try:	
                ##could be two values then always bigger like in gnome project
                val = float(val)	
    			
                if(val>=0.5 and int(cs[0])==0) or (val<0.5 and int(cs[0])==1):
                    tmp=[e[1],e[2],e[3],e[4],val]
                    gnome.append(tmp)
                    individuals= individuals + 1
            except ValueError:
                if(int(cs[0])==1):
                    individuals= individuals + 1
                    tmp=[e[1],e[2],e[3],e[4],val]
                    gnome.append(tmp)
        gnomes.append(gnome)
    print("Found " + str(individuals) + " individuals!")
    return gnomes

def createTarWithGnomes(gnomes,columns,c):
    #print(gnomes)
    tarfile_name = "/tmp/file.tar.gz"
    tar = tarfile.open(name = tarfile_name, mode = "w:gz")
    for i in range(len(gnomes)):
        for j in range(len(gnomes[i])):
            gnomes[i][j] = "\t".join(str(e) for e in gnomes[i][j])
        gnomes[i]="\n".join(gnomes[i])

        string = io.BytesIO(gnomes[i].encode())
        info = tarfile.TarInfo(name="chr"+str(c)+"."+columns[i])
        info.size = len(string.getbuffer())
        tar.addfile(tarinfo=info, fileobj=string)
    tar.close()
    return tarfile_name 

def readGzip(key_datafile, inspector):
    #Read data
    tmp_file = "/tmp/file.vcf.gz"
    
    inspector.addTimeStamp("Time Stamp before Download 2")
    
    pyStorage.copy(key_datafile, tmp_file)
    
    inspector.addTimeStamp("Time Stamp after Download 2")
    
    
    with gzip.open(tmp_file, 'rb') as f:
        file_content = f.read().decode("utf-8") 
    
    return file_content.split('\n')
    
    
def readColumns(key_columnsfile, inspector):
    ##Read columns
    tmp_file = "/tmp/file.txt"
    
    inspector.addTimeStamp("Time Stamp before Download 1")
    pyStorage.copy(key_columnsfile, tmp_file)
    inspector.addTimeStamp("Time Stamp after Download 1")
    
    object = open(tmp_file, "r").read()
    bytes(object, "utf-8")

    return object.rstrip().split("\t")
