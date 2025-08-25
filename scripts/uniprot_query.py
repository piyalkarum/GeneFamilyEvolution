from bioservices import UniProt
import csv

# Initialize the UniProt service
uniprot = UniProt()

# Function to fetch annotations for a single gene ID
def fetch_gene_details(gene_id):
    try:
        # Retrieve results in tabular format
        result = uniprot.search(gene_id)
        
        if result:
            # Parse the TSV data into a list of dictionaries
            lines = result.strip().split("\n")
            headers = lines[0].split("\t")
            data = [dict(zip(headers, line.split("\t"))) for line in lines[1:]]
            
            # Add gene_id to each row
            for row in data:
                row["gene_id"] = gene_id
            return data
        else:
            print(f"No data found for {gene_id}")
            return []
    except Exception as e:
        print(f"Error fetching data for {gene_id}: {e}")
        return []

# Read gene IDs from a file
input_file = "gene_ids_to_search_functions.txt"  # Replace with file name
with open(input_file, "r") as file:
    gene_ids = [line.strip() for line in file]

# Fetch details for all gene IDs
results = []
for gene_id in gene_ids:
    print(f"Fetching data for gene ID: {gene_id}")
    data = fetch_gene_details(gene_id)
    results.extend(data)

# Save results to a CSV file
output_file = "gene_function_annotations.csv"
if results:
    # Extract headers from the first result
    headers = list(results[0].keys())  # Include all keys, including "gene_id"
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        writer.writerows(results)
    print(f"Annotations saved to {output_file}")
else:
    print("No results retrieved.")