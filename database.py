import sqlalchemy

common_connection = None
preset_database_host = "localhost"
preset_database_port = 5000

def create_database_connection(host=preset_database_host, port=preset_database_port):
    global common_connection
    if common_connection is None: 
        engine = sqlalchemy.create_engine(
            "postgresql://postgres:postgres@%s:%d/bioinfo" % (host, port))
        common_connection = engine.connect()

    return common_connection

def close_database_connection():
    global common_connection
    if common_connection is not None:
        common_connection.close()
        common_connection = None
