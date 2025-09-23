#!/usr/bin/env python3
import os
import glob
import subprocess
import sys

def find_rebuilt_chunks(rebuilt_dir):
    """Find all rebuilt chunk files"""
    return sorted(glob.glob(os.path.join(rebuilt_dir, "*_rebuilt.h5ad")))

def submit_aucell_job(rebuilt_file, regulons_csv, output_dir):
    """Submit AUCell job for a rebuilt chunk"""
    base_name = os.path.basename(rebuilt_file).replace("_rebuilt", "")
    aucell_output = os.path.join(output_dir, f"{base_name}_aucell.h5ad")
    
    cmd = [
        "sbatch",
        "60_aucell_patched.sh",
        rebuilt_file,
        regulons_csv,
        aucell_output,
    ]
    return cmd, aucell_output

def main():
    if len(sys.argv) != 2:
        print("Usage: python submit_aucell_jobs.py <rebuilt_chunks_directory>")
        sys.exit(1)
    
    rebuilt_dir = sys.argv[1]
    regulons_csv = "/scratch/easmit31/GRN_copy/scenic/h5ad_files/astros_allLv_regulons.csv"
    
    aucell_dir = rebuilt_dir.replace("_rebuilt", "_aucell")
    os.makedirs(aucell_dir, exist_ok=True)
    
    chunks = find_rebuilt_chunks(rebuilt_dir)
    print(f"Found {len(chunks)} rebuilt chunks to process")
    
    job_ids = []
    
    for i, rebuilt_file in enumerate(chunks, 1):
        print(f"\n=== Submitting job {i}/{len(chunks)}: {os.path.basename(rebuilt_file)} ===")
        
        aucell_cmd, aucell_output = submit_aucell_job(rebuilt_file, regulons_csv, aucell_dir)
        print("Command:", " ".join(aucell_cmd))
        
        result = subprocess.run(aucell_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            print(f"✅ Submitted job {job_id}")
        else:
            print(f"❌ ERROR submitting job for {rebuilt_file}")
            print(result.stderr)
    
    print("\n=== SUBMISSION COMPLETE ===")
    print(f"Submitted {len(job_ids)} AUCell jobs")
    for jid in job_ids:
        print(f"  {jid}")
    
    print("\nMonitor with: squeue -u $USER")
    print(f"Results will be in: {aucell_dir}")

if __name__ == "__main__":
    main()
