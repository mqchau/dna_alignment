import sqlalchemy
import os

common_connection = None
preset_database_host = os.environ["RDS_HOST"]
preset_database_port = int(os.environ["RDS_PORT"])
preset_database_name = os.environ["RDS_DATABASE"]

def create_database_connection(host=preset_database_host, port=preset_database_port, database=preset_database_name):
    global common_connection
    if common_connection is None: 
        engine = sqlalchemy.create_engine(
            "postgresql://postgres:postgres@%s:%d/%s" % (host, port, database))
        common_connection = engine.connect()

    return common_connection

def close_database_connection():
    global common_connection
    if common_connection is not None:
        common_connection.close()
        common_connection = None
