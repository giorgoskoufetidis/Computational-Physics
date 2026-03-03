import subprocess
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
app_path = os.path.join(HERE, "star_app.py")

subprocess.run([
    sys.executable,
    "-m", "streamlit",
    "run", 
    app_path,
    "--server.headless=true",
    "--browser.gatherUsageStats=false"
])
