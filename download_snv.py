cohort = 'PAAD'
import requests
import json
import os
import gzip
import time
from urllib.request import urlretrieve
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import RequestException

def create_session_with_retry():
    """Create a requests session with retry strategy"""
    session = requests.Session()
    retries = Retry(
        total=5,  # number of retries
        backoff_factor=1,  # wait 1, 2, 4, 8, 16 seconds between retries
        status_forcelist=[500, 502, 503, 504, 104],  # HTTP status codes to retry on
        allowed_methods=["HEAD", "GET", "OPTIONS"]  # HTTP methods to retry
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

def make_request_with_retry(url, params=None, headers=None, max_retries=5):
    """Make HTTP request with retry logic"""
    session = create_session_with_retry()
    
    for attempt in range(max_retries):
        try:
            if params:
                response = session.get(url, params=params, headers=headers)
            else:
                response = session.get(url, headers=headers)
            response.raise_for_status()
            return response
        except RequestException as e:
            if attempt == max_retries - 1:  # Last attempt
                raise e
            wait_time = (2 ** attempt) + 1  # Exponential backoff
            print(f"Request failed, retrying in {wait_time} seconds... (Attempt {attempt + 1}/{max_retries})")
            time.sleep(wait_time)

def get_files():
    files_endpt = "https://api.gdc.cancer.gov/files"
    filters = {
        "op": "and",
        "content":[
            {
                "op": "=",
                "content":{
                    "field": "cases.project.project_id",
                    "value": "TCGA-"+cohort
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "data_type",
                    "value": "Masked Somatic Mutation"
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "file_name",
                    "value": "*masked.maf.gz*"
                }
            }
        ]
    }
   
    params = {
        "filters": json.dumps(filters),
        "fields":"file_id,file_name,cases.submitter_id",
        "format": "JSON",
        "size": "1000"
    }
   
    response = make_request_with_retry(files_endpt, params=params)
    file_data = json.loads(response.content.decode("utf-8"))
   
    return file_data["data"]["hits"]

def download_and_process_files(file_hits, output_dir):
    data_endpt = "https://api.gdc.cancer.gov/data/"
   
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
   
    for hit in file_hits:
        file_id = hit["file_id"]
        file_name = hit["file_name"]
        case_id = hit["cases"][0]["submitter_id"]
       
        temp_gz_file = os.path.join(output_dir, f"{file_name}")
        final_output_file = os.path.join(output_dir, f"{case_id}.maf")
        if os.path.exists(final_output_file):
            continue
       
        time.sleep(0.1)
        try:
            response = make_request_with_retry(
                f"{data_endpt}{file_id}", 
                headers={"Content-Type": "application/json"}
            )
            
            with open(temp_gz_file, 'wb') as f:
                f.write(response.content)
           
            with gzip.open(temp_gz_file, 'rb') as gz_file:
                with open(final_output_file, 'wb') as output_file:
                    output_file.write(gz_file.read())
           
            os.remove(temp_gz_file)
            print(f"Processed {case_id}.maf")
            
        except Exception as e:
            print(f"Error processing {case_id}: {str(e)}")
            if os.path.exists(temp_gz_file):
                os.remove(temp_gz_file)
            continue

def main():
    print("Exploring files...")
    try:
        file_hits = get_files()
       
        if not file_hits:
            print("No files found matching criteria")
            return
       
        print(f"Found {len(file_hits)} matching files")
       
        print("Starting downloads...")
        download_and_process_files(file_hits, "data_snv_"+cohort)
        print("Download and processing complete!")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()