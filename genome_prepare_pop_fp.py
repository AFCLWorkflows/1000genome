from Inspector import *

def lambda_handler(event, context):

    inspector = Inspector()
    inspector.inspectAll()
    
    inspector.addAttribute("function_name", "genome_pop_nv")
    inspector.addTimeStamp("Time Stamp at start")
    
    failure_prob = event['failure_prob']
    output_buckets = event["output_buckets"]
    key_columnsfile = event['key_columnsfile']
    siftfile = event['sifted']
    merged_result = event['merged_result']
    
    
    inspector.addTimeStamp("Time Stamp at end")
    
    inspector.addAttribute("POP",  [event['AFR'],event['ALL'],event['AMR'],event['EAS'],event['EUR'],event['GBR'],event['SAS']])
    inspector.addAttribute("failure_prob", failure_prob)
    inspector.addAttribute("output_buckets", output_buckets)
    inspector.addAttribute("sifted", siftfile)
    inspector.addAttribute("key_columnsfile", key_columnsfile)
    inspector.addAttribute("merged_result", merged_result)
    inspector.addAttribute("aws_access_key_id", event['aws_access_key_id'])
    inspector.addAttribute("aws_secret_key", event['aws_secret_key'])
    inspector.addAttribute("aws_session_token", event['aws_session_token'])
    inspector.addAttribute("gcp_client_email", event['gcp_client_email'])
    inspector.addAttribute("gcp_private_key", event['gcp_private_key'])
    inspector.addAttribute("gcp_project_id", event['gcp_project_id'])
    inspector.inspectAllDeltas()
    return inspector.finish() 