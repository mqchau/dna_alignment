import sqlalchemy
import os

common_connection = {}
preset_database_host = os.environ["RDS_HOST"]
preset_database_port = int(os.environ["RDS_PORT"])
preset_database_name = os.environ["RDS_DATABASE"]

def create_database_connection(host=preset_database_host, port=preset_database_port, database=preset_database_name):
    global common_connection
    if database not in common_connection: 
        engine = sqlalchemy.create_engine(
            "postgresql://postgres:postgres@%s:%d/%s" % (host, port, database))
        common_connection[database] = engine.connect()
        

    return common_connection[database]

def close_database_connection(database=preset_database_name):
    global common_connection
    if database in common_connection: 
        common_connection[database].close()
        del common_connection[database]
