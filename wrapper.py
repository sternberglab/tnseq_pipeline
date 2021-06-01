import botocore
import boto3
import os
import json
import shutil
from pathlib import Path

from main import main

def convert_dynamo_item_to_json(dyn_item):
	new_item = dict()
	for k,v in dyn_item.items():
		if 'S' in v:
			new_item[k] = v['S']
		elif 'BOOL' in v:
			new_item[k] = not not v['BOOL']
		elif 'N' in v:
			new_item[k] = int(v['N'])
	return new_item

def wrapper():
	with open('./sqs_message.json', 'r') as info_file:
		info = json.load(info_file)
	sample = info['Sample']
	analysisId = info['analysisId']

	dynamo = boto3.client('dynamodb')
	paginator = dynamo.get_paginator('query')
	for page in paginator.paginate(
		TableName='Illumina_dataset',
		IndexName='appsession_id',
		KeyConditionExpression="appsession_id=:analysisId",
		ExpressionAttributeValues={':analysisId': {'S': analysisId}, ':sample': {'S': sample}},
		ExpressionAttributeNames={'#name': 'name'},
		FilterExpression='#name=:sample',
		Limit=1000
	):
		if len(page['Items']):
			item = page['Items'][0]

	s3 = boto3.client('s3')
 
	metadata = main(isCloud=True)
	dynamo.update_item(
		TableName='Illumina_dataset',
		Key={
			'id': item['id']
		},
		UpdateExpression="SET pipeline_info=:pipeline_info",
		ExpressionAttributeValues={
			':pipeline_info': {'S': json.dumps(metadata)},
		}
	)
	
	s3prefix = f"ngs_pipeline_outputs/{analysisId}/{sample}"
	genome_locations_filepath = f"./outputs/samples/{sample}_genome_read_locations.csv"
	s3.upload_file(genome_locations_filepath, 'sternberg-sequencing-data', f"{s3prefix}/genome_read_locations.csv")
	plots_filestem = f"./outputs/plots/{sample}"
	s3.upload_file(f"{plots_filestem}_genome_hist_normalized.svg", 'sternberg-sequencing-data', f"{s3prefix}/plots/genome_normalized.svg")
	s3.upload_file(f"{plots_filestem}_genome_hist_raw.svg", 'sternberg-sequencing-data', f"{s3prefix}/plots/genome_raw.svg")
	s3.upload_file(f"{plots_filestem}_genome_hist_zoomed.svg", 'sternberg-sequencing-data', f"{s3prefix}/plots/genome_zoomed.svg")
	s3.upload_file(f"{plots_filestem}_trans_dist_hist_overlap_v2.svg", 'sternberg-sequencing-data', f"{s3prefix}/plots/trans_dist_overlap.svg")
	s3.upload_file(f"{plots_filestem}_trans_dist_hist_no_overlap_v2.svg", 'sternberg-sequencing-data', f"{s3prefix}/plots/trans_dist_no_overlap.svg")
	
	# remove the tmp directory with raw files and intermediates incase this container needs to process many samples
	shutil.rmtree(os.path.join(Path(__file__).parent.absolute(), 'tmp'))
wrapper()