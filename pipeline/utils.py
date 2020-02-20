from pathlib import Path
import os

intermediates_dir = os.path.join(Path(__file__).cwd(), 'intermediates')
os.makedirs(intermediates_dir, exist_ok=True)

def inter_path(filename):
	return os.path.join(intermediates_dir, filename)
