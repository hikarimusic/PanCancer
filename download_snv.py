cohort = 'PAAD'

import requests
import json
import os
import gzip
import time
from urllib.request import urlretrieve

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
    
    response = requests.get(files_endpt, params=params)
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
        
        time.sleep(0.1)
        response = requests.get(f"{data_endpt}{file_id}", headers={"Content-Type": "application/json"})
        with open(temp_gz_file, 'wb') as f:
            f.write(response.content)
        
        with gzip.open(temp_gz_file, 'rb') as gz_file:
            with open(final_output_file, 'wb') as output_file:
                output_file.write(gz_file.read())
        
        os.remove(temp_gz_file)
        print(f"Processed {case_id}.maf")

def main():
    print("Exploring files...")
    file_hits = get_files()
    
    if not file_hits:
        print("No files found matching criteria")
        return
    
    print(f"Found {len(file_hits)} matching files")
    
    print("Starting downloads...")
    download_and_process_files(file_hits, "data_snv_"+cohort)
    print("Download and processing complete!")

if __name__ == "__main__":
    main()