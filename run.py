"""
run.py

Automates the setup and execution of the crop productivity analysis pipeline.
Assumes you are in the root of the crop_productivity_project.
"""

import os
import subprocess
import sys

def run_shell(cmd):
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"Error: Command failed with exit code {result.returncode}")
        sys.exit(result.returncode)

def main():
    print("\n🔧 STEP 1: Create virtual environment")
    run_shell("python -m venv venv")

    print("\n✅ STEP 2: Activate virtual environment and install dependencies")
    activate_cmd = "source venv/bin/activate" if os.name != 'nt' else "venv\\Scripts\\activate"
    pip_install_cmd = f"{activate_cmd} && pip install -r requirements.txt && pip install -e ."

    print("\n📦 Installing packages...")
    run_shell(f"/bin/bash -c '{pip_install_cmd}'")

    print("\n🚀 STEP 3: Run the CLI pipeline")
    run_shell(f"/bin/bash -c '{activate_cmd} && cropproductivity --params cropProductivity.par --output GEE_ETH'")

if __name__ == "__main__":
    main()
