from pathlib import Path
import pandas as pd

manifest_data = []

# Loop over all subfolders in raw_data
for sample_folder in raw_dir.iterdir():
    if sample_folder.is_dir():
        sample_id = sample_folder.name
        files = list(sample_folder.glob("*"))  # all files in that folder
        
        if sequence_type == "NGS" and region == "V3V4":  # paired-end
            for f in files:
                direction = "forward" if "_R1" in f.stem else "reverse"
                manifest_data.append({
                    "sampleid": sample_id,
                    "absolute-filepath": f.resolve(),
                    "direction": direction
                })
        elif sequence_type == 'TGS' and region == "Full" : # single-end full-length
            for f in files:
                manifest_data.append({
                    "sampleid": sample_id,
                    "absolute-filepath": f.resolve()
                })
        else:
            raise ValueError("Unsupported sequence type or region configuration.")

# Write manifest
manifest_file = study_dir / "manifest.tsv"
pd.DataFrame(manifest_data).to_csv(manifest_file, index=False), sep='\t'
print(f"Manifest written to {manifest_file}")