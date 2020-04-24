from appJar import gui
import csv
from pathlib import Path
from pipeline.main import main

def get_genomes():
	genomes = []
	genomes_csv = Path('../resources/genomes/genomes.csv')
	if not genomes_csv.exists():
		raise
	with open(genomes_csv, 'r', encoding='utf-8-sig') as file:
		reader = csv.DictReader(file)
		genomes = [row for row in reader]
	return genomes

def get_systems():
	systems = []
	systems_csv = Path('../resources/systems/systems.csv')
	if not systems_csv.exists():
		raise
	with open(systems_csv, 'r', encoding='utf-8-sig') as file:
		reader = csv.DictReader(file)
		systems = [row for row in reader]
	return systems

#genomes = get_genomes()
#systems = get_systems()

def run():
	values = app.getAllInputs()
	if "input" not in values:
		app.setLabel("inputError", "Must provide a file")
		app.setLabelFg("inputError", "red")
		return
	else:
		app.setLabel("inputError", "")
	#genome = next((item for item in genomes if item["name"] == values["genome"]))
	#system = next((item for item in systems if item["name"] == values["system"]))
	info_file = values["input"]
	main(info_file)
	print("finished it all, back to gui")

app = gui("Sternberg lab tools")

app.setSize(500,500)

app.addLabel("title", "Illumina pipeline")
#app.addLabelOptionBox("genome", [genome["name"] for genome in get_genomes()])
#app.addLabelOptionBox("system", [system["name"] for system in get_systems()])
app.addLabel("inputError", "")
app.addFileEntry("input")
app.addButton("Run", run)

app.go()