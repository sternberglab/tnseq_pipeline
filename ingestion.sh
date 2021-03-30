#!/bin/bash

echo "ECS Illumina pipeline running"

git pull
$(aws secretsmanager get-secret-value --secret-id NGS_PIPELINE_UNPROCESSED_SQS_URL |  jq -r '"export NGS_PIPELINE_UNPROCESSED_SQS_URL=" + .SecretString ')

while [ /bin/true ]; do

	msg=$( \
      aws sqs receive-message \
          --queue-url $NGS_PIPELINE_UNPROCESSED_SQS_URL \
          --wait-time-seconds 20 \
          --output text \
          --query Messages[0].[Body,ReceiptHandle] \
          --visibility-timeout 300
    )

    if [ -z "${msg}" -o "${msg}" = "None" ]; then
    	echo "$(date) Processing complete. Stopping task."
    	exit
    else
		echo "$(date) SQS Message: ${msg}"
		sqs_message=$(echo "${msg}" | cut -f1 --)
		echo "${sqs_message}" > ./sqs_message.json

		receipt_handle=$(echo "${msg}" | cut -f2 --)
		$(python3 ./test.py)
		CMD_EXIT=$?

		if [ $CMD_EXIT -eq 0 ]; then
			echo "SUCCESS"
		else
			echo "FAILURE"
		fi

		rm ./sqs_message.json

		aws sqs delete-message \
			--queue-url $NGS_PIPELINE_UNPROCESSED_SQS_URL \
			--receipt-handle ${receipt_handle}
	fi
done
exit
