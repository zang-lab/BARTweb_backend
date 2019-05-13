#!/bin/bash

# Requires the aws-cli tools installed
# Follow these steps: https://gist.github.com/nmagee/0ce2cd44a1c4664783a68a903a4d6264
# 
# Requires the jq JSON parser
# Follow these steps: https://stedolan.github.io/jq/

BART_dir="/sfs/bart"
LOG_dir=$BART_dir/log

# Check the queue for message count.
QCOUNT=`/home/wm9tr/miniconda3/bin/aws sqs get-queue-attributes \
  --queue-url "https://sqs.us-east-1.amazonaws.com/474683445819/bart-web-processing" \
  --attribute-names "ApproximateNumberOfMessages" | \
  /home/wm9tr/bin/jq -r .Attributes.ApproximateNumberOfMessages`

if [ "$QCOUNT" -gt 0 ]; then
  current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
  echo $current_date_time >> $LOG_dir/aws_queue.log
  echo $QCOUNT >> $LOG_dir/aws_queue.log
  # Get a message from the queue if there is more than 0 message
  RAW=`/home/wm9tr/miniconda3/bin/aws sqs receive-message \
    --message-attribute-names "submissionkey" \
    --max-number-of-messages 1 \
    --queue-url "https://sqs.us-east-1.amazonaws.com/474683445819/bart-web-processing" \
    --wait-time-seconds 20`;

  # The $DIR variable is the name of your submissionkey
  DIR=`echo $RAW | /home/wm9tr/bin/jq -r .Messages[0].MessageAttributes.submissionkey.StringValue`;

  # execute the trigger.py
  /home/wm9tr/miniconda3/bin/python $BART_dir/trigger.py $DIR >> $LOG_dir/aws_queue.log 2>&1;
  
  # Finally if you have gotten the message and worked with it, you must delete it or it will remain in the queue
  RECEIPTHANDLE=`echo $RAW | /home/wm9tr/bin/jq -r .Messages[0].ReceiptHandle`;
  echo $RECEIPTHANDLE >> $LOG_dir/aws_queue.log;
  /home/wm9tr/miniconda3/bin/aws sqs delete-message \
    --queue-url "https://sqs.us-east-1.amazonaws.com/474683445819/bart-web-processing" \
    --receipt-handle "$RECEIPTHANDLE";  

else

  # do nothing - if this script is on a cron job and runs every 1 minute, this will happen most of the time
  # or log the check if you want
  echo "There is nothing to do.";
fi
