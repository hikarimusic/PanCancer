import pandas as pd
import os
import re

def combine_clinical_cohorts(directory_path):
    """
    Combines clinical patient data from multiple cancer cohorts into a single DataFrame.
    
    Parameters:
    directory_path (str): Path to directory containing the clinical patient files
    
    Returns:
    pandas.DataFrame: Combined data from all cohorts with cohort information
    """
    # Initialize empty list to store DataFrames
    dfs = []
    
    # Pattern to match clinical patient files and extract cancer type
    pattern = r'clinical_patient_([a-zA-Z]+)\.txt$'
    
    # Iterate through files in directory
    for filename in os.listdir(directory_path):
        if 'clinical_patient' in filename and filename.endswith('.txt'):
            file_path = os.path.join(directory_path, filename)
            
            # Extract cancer type from filename
            match = re.search(pattern, filename)
            if match:
                cancer_type = match.group(1).upper()
                
                try:
                    # Read the file
                    # Using sep='\t' since the files appear to be tab-delimited
                    # Read headers first (first row)
                    headers = pd.read_csv(file_path, sep='\t', nrows=1)
                    
                    # Read the actual data (skip the first 3 rows - header + 2 metadata rows)
                    df = pd.read_csv(file_path, sep='\t', skiprows=[1, 2], names=headers.columns)
                    
                    # Add cancer type column
                    df['cohort'] = cancer_type
                    
                    # Add source filename
                    df['source_file'] = filename
                    
                    dfs.append(df)
                    print(f"Successfully processed {filename}")
                    
                except Exception as e:
                    print(f"Error processing {filename}: {str(e)}")
    
    if not dfs:
        raise ValueError("No valid files were processed")
    
    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Basic cleaning
    # Remove any completely empty columns
    combined_df = combined_df.dropna(axis=1, how='all')
    
    # Remove any duplicate rows
    combined_df = combined_df.drop_duplicates()
    
    # Print summary statistics
    print("\nSummary:")
    print(f"Total number of patients: {len(combined_df)}")
    print("\nPatients per cohort:")
    print(combined_df['cohort'].value_counts())
    print("\nColumns in combined dataset:", len(combined_df.columns))
    
    return combined_df

# Example usage:
if __name__ == "__main__":
    # Replace with your directory path
    directory_path = "./data_clinical"
    
    try:
        combined_data = combine_clinical_cohorts(directory_path)
        
        # Save combined data
        output_file = "combined_clinical_cohorts.csv"
        combined_data.to_csv(output_file, index=False)
        print(f"\nCombined data saved to {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")