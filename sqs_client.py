import boto3
import os

sqs = boto3.resource('sqs')

all_queues = {
    "reference_hash" : sqs.get_queue_by_name(QueueName='cm124_reference_hash_queue'),
    "align_base": sqs.get_queue_by_name(QueueName='cm124_align_base_queue'),
    "get_mutation": sqs.get_queue_by_name(QueueName='cm124_get_mutation_queue')
}

def send_message(queue_name, message):
    global all_queues
    all_queues[queue_name].send_message(MessageBody=message)

def receive_message(queue_name):
    global all_queues
    return all_queues[queue_name].receive_messages(MaxNumberOfMessages=1)[0]


if __name__ == "__main__":
    send_message('reference_hash', '0,1000')
    new_message=receive_message('reference_hash')
    print new_message.body
    new_message.delete()
