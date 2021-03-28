#!/bin/bash

echo "ECS Illumina pipeline running"

$(aws ssm get-parameters --with-decryption --names NGS_PIPELINE_UNPROCESSED_SQS_URL |  jq -r '.Parameters| .[] | "export " + .Name + "=" + .Value ')

echo "got thing $(NGS_PIPELINE_UNPROCESSED_SQS_URL)"

while [ /bin/true ]; do
	$(git pull)

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
		echo "${sqs_message}" > Illumina-pipeline/sqs_message.json

		receipt_handle=$(echo "${msg}" | cut -f2 --)
		$(python3 Illumina-pipeline/test.py)
		CMD_EXIT=$?

		if [ $CMD_EXIT -eq 0 ]; then
			echo "SUCCESS"
		else
			echo "FAILURE"
		fi

		rm Illumina-pipeline/sqs_message.json

		aws sqs delete-message \
			--queue-url $NGS_PIPELINE_UNPROCESSED_SQS_URL \
			--receipt-handle ${receipt_handle}
	fi
done
exit
