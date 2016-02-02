import ipdb
import boto3
import os

sqs = boto3.resource('sqs')

all_queues = {
    "reference_hash" : sqs.get_queue_by_name(QueueName='cm124_reference_hash_queue'),
    "align_base": sqs.get_queue_by_name(QueueName='cm124_align_base_queue'),
    "get_mutation": sqs.get_queue_by_name(QueueName='cm124_get_mutation_queue'),
    "start_job": sqs.get_queue_by_name(QueueName='cm124_start_job_queue'),
}

def send_message(queue_name, message):
    global all_queues
    all_queues[queue_name].send_message(MessageBody=message)

def receive_message(queue_name):
    global all_queues
    messages = all_queues[queue_name].receive_messages(MaxNumberOfMessages=1)
    if len(messages) > 0:
        return messages[0]
    return None

def get_queue_remaining_messages(queue_name):
    global all_queues
    all_queues[queue_name].load()
    # ipdb.set_trace()
    return int(all_queues[queue_name].attributes['ApproximateNumberOfMessages']) + int(all_queues[queue_name].attributes['ApproximateNumberOfMessagesNotVisible'])

if __name__ == "__main__":
    # send_message('reference_hash', '0,1000')
    new_message=receive_message('reference_hash')
    if new_message is not None:
        print new_message.body
        print get_queue_remaining_messages('reference_hash')
        new_message.delete()
