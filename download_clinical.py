import requests
import json
import os
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
                    "value": "*TCGA*"
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "data_type",
                    "value": "Clinical Supplement"
                }
            },
            {
                "op": "=",
                "content":{
                    "field": "file_name",
                    "value": "*clinical_patient*"
                }
            }
        ]
    }

    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name",
        "format": "JSON",
        "size": "100"
    }
    
    response = requests.get(files_endpt, params=params)
    file_data = json.loads(response.content.decode("utf-8"))
    
    return file_data["data"]["hits"]


def download_files(file_hits, output_dir):
    data_endpt = "https://api.gdc.cancer.gov/data/"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for hit in file_hits:
        file_id = hit["file_id"]
        file_name = hit["file_name"]
        output_file = os.path.join(output_dir, f"{file_name}")
        
        print(f"Downloading {file_name}...")
        response = requests.get(f"{data_endpt}{file_id}", headers={"Content-Type": "application/json"})
        with open(output_file, 'wb') as f:
            f.write(response.content)


def main():
    print("Exploring files ...")
    file_hits = get_files()
    # for file in file_hits:
    #     print(file)
    
    if not file_hits:
        print("No files found matching")
        return
    print(f"Found {len(file_hits)} matching files")
    
    print("Starting downloads...")
    download_files(file_hits, "data_clinical")
    print("Download complete!")

if __name__ == "__main__":
    main()



# file_id
# file_name
# cases.submitter_id
# cases.case_id
# data_category
# data_type
# cases.samples.tumor_descriptor
# cases.samples.tissue_type
# cases.samples.sample_type
# cases.samples.submitter_id
# cases.samples.sample_id
# analysis.workflow_type
# cases.project.project_id
# cases.samples.portions.analytes.aliquots.aliquot_id
# cases.samples.portions.analytes.aliquots.submitter_id