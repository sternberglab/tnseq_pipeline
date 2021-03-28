print("THIS IS A TEST")
import botocore
import boto3

def main():
	#with open('sqs-message.json', 'r') as msg:
	s3 = boto3.client('s3')
	s3.upload_file('./sqs-message.json', 'sternberg-sequencing-data', 'test/test.py')
 
main()