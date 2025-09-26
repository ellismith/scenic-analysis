#!/scratch/easmit31/conda_envs/pyscenic_v2/bin/python
import os
import glob
import subprocess
import sys

def find_rebuilt_chunks(rebuilt_dir):
    """Find all rebuilt chunk files"""
    return sorted(glob.glob(os.path.join(rebuilt_dir, "*_rebuilt.h5ad")))

def submit_aucell_job(rebuilt_file, regulons_csv, output_dir):
    """Submit AUCell job for a rebuilt chunk"""
    base_name = os.path.basename(rebuilt_file).replace("_rebuilt.h5ad", "_aucell.h5ad")
    aucell_output = os.path.join(output_dir, base_name)
    
    cmd = [
        "sbatch",
        "60_aucell_patched.sh",
        rebuilt_file,
        regulons_csv,
        aucell_output,
    ]
    return cmd, aucell_output

def main():
    if len(sys.argv) != 3:
        print("Usage: python submit_aucell_jobs.py <rebuilt_chunks_directory> <regulons_csv>")
        print("\nExample:")
        print("  python submit_aucell_jobs.py /scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_adults_rebuilt /scratch/easmit31/GRN_copy/scenic/h5ad_files/gaba_allLv_regulons.csv")
        sys.exit(1)
    
    rebuilt_dir = sys.argv[1]
    regulons_csv = sys.argv[2]
    
    # Verify inputs exist
    if not os.path.isdir(rebuilt_dir):
        print(f"ERROR: Rebuilt directory not found: {rebuilt_dir}")
        sys.exit(1)
    
    if not os.path.isfile(regulons_csv):
        print(f"ERROR: Regulons CSV not found: {regulons_csv}")
        sys.exit(1)
    
    # Create output directory
    aucell_dir = rebuilt_dir.replace("_rebuilt", "_aucell")
    os.makedirs(aucell_dir, exist_ok=True)
    
    chunks = find_rebuilt_chunks(rebuilt_dir)
    
    if not chunks:
        print(f"ERROR: No rebuilt chunks found in {rebuilt_dir}")
        sys.exit(1)
    
    print(f"Found {len(chunks)} rebuilt chunks to process")
    print(f"Regulons CSV: {regulons_csv}")
    print(f"Output directory: {aucell_dir}")
    
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
