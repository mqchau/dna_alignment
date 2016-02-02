import sqs_client

if __name__ == "__main__":
    sqs_client.send_message('start_job', 'reference_hash,practice')
