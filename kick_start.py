import sqs_client

if __name__ == "__main__":
    # sqs_client.send_message('start_job', 'align_base,practice2')
    sqs_client.send_message('start_job', 'get_mutation,practice2')
