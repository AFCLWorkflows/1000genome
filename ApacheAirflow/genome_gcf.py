from datetime import timedelta
import json

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG
from airflow.models import Variable
# Operators; we need this to operate!
from airflow.operators.bash_operator import BashOperator
from airflow.operators.http_operator import SimpleHttpOperator
from airflow.operators.python_operator import PythonOperator
from airflow.operators.dummy_operator import DummyOperator
from airflow.utils.dates import days_ago
from airflow.operators.subdag_operator import SubDagOperator
#from airflow.models import get_fernet
# These args will get passed on to each operator
# You can override them on a per-task basis during operator initialization
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'email': ['airflow@example.com'],
    'email_on_failure': False,
    'email_on_retry': False,
    'start_date': days_ago(2),
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}


"""
Preparation:

1. Add gcf_conn to Admin->Connections section of the Airflow UI.
   Consists of the HTTP Trigger URL without the endpoint (https://your-region-your-project-name.cloudfunctions.net)

"""


DAG_NAME='genome_gcf'

dag = DAG(
    DAG_NAME,
    default_args=default_args,
    description='Genome Workflow to determine mutation overlap and frequency',
    schedule_interval=None,

)


start = DummyOperator(
	task_id = 'start',
	dag = dag

)


#Note: Input can also be read from input json and stored to Airflow Variables and/or XCom. Outputs can be passed on to other functions via XCom.
key = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz"
key_input = "50_individuals.vcf.gz"
k = 5
x = [0,1,2,3,4]



sifting = SimpleHttpOperator(
    task_id='sifting',
    method='POST',
    http_conn_id='gcf_conn',
    endpoint='genome_sifting',
    headers={"Accept":"application/json", "Content-Type": "application/json"},
    data= json.dumps({'key': key}),
    dag=dag
)


def prepare_individuals():
	
	finalstring = "{\"key_datafiles\": ["
	for i in range(0, k):
		string = "\"chr22n-" + str(int((i*10000/k)+1)) + "-" + str(int((i+1)*(10000/k)+1)) + ".tar.gz\""
		print(string)
		finalstring = finalstring + string 
		if i != k-1:
			finalstring = finalstring + ", "
	
	finalstring = finalstring + "] }"	
	print(finalstring)
	return finalstring	
		

prepare_merge_input = PythonOperator(
	task_id='prepare_merge_input',
	python_callable=prepare_individuals,
	dag=dag

)




prepare_pop =  DummyOperator(
	task_id='prepare_pop',
	dag=dag
)



for i in range(int(k)):
	individual = SimpleHttpOperator(
		task_id='individual_' + str(i),
		method='POST',
		http_conn_id='gcf_conn',
		endpoint='genome_individual',
                headers={"Accept":"application/json", "Content-Type": "application/json"},
		data=json.dumps({'key_input': key_input, 'key': key, 'k': int(k), 'x': i}),
		dag=dag
	)

	start.set_downstream(individual)
	individual.set_downstream(prepare_merge_input)



individuals_merge = SimpleHttpOperator(
	task_id='individuals_merge',
	method='POST',
	http_conn_id='gcf_conn',
	endpoint='genome_individuals_merge',
	headers={"Accept":"application/json", "Content-Type": "application/json"},
	data = json.dumps({"key_datafiles": ["chr22n-1-2001.tar.gz", "chr22n-2001-4001.tar.gz", "chr22n-4001-6001.tar.gz", "chr22n-6001-8001.tar.gz", "chr22n-8001-10001.tar.gz"]}),
	dag = dag
	
	)




populations=["AFR", "ALL", "AMR", "EAS", "EUR", "GBR", "SAS"]


for i in range(len(populations)):
	frequency = SimpleHttpOperator(
		task_id='frequency_' + str(i),
		method='POST',
		http_conn_id='gcf_conn',
		endpoint='genome_frequency',
        headers={"Accept":"application/json", "Content-Type": "application/json"},
		data=json.dumps({'POP': populations[i]}),
		dag=dag
	)
	
	prepare_pop.set_downstream(frequency)
	frequency.set_upstream(individuals_merge)



for i in range(len(populations)):
	mutual_overlap = SimpleHttpOperator(
		task_id='mutual_overlap_' + str(i),
		method='POST',
		http_conn_id='gcf_conn',
		endpoint='genome_mutual_overlap',
        headers={"Accept":"application/json", "Content-Type": "application/json"},
		data=json.dumps({'POP': populations[i]}),
		dag=dag
	)
	
	prepare_pop.set_downstream(mutual_overlap)
	mutual_overlap.set_upstream(individuals_merge)


start>>sifting

prepare_merge_input>>prepare_pop

prepare_merge_input>>individuals_merge



