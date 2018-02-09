from helioviewer.db  import *
import os
db, cursor = setup_database_schema(os.environ["MYSQL_ROOT_USER"], os.environ["MYSQL_ROOT_PASSWORD"],
            "swhv-mysql", os.environ["HVDBNAME"], os.environ["MYSQL_ROOT_USER"], os.environ["MYSQL_ROOT_PASSWORD"], True)
