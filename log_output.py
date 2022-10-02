import logging
from os import mkdir, path, remove


""" 
The results or Logger can be summarized in Log/ folder
Each Logging can be turned off by 'lg.setlevel(logging.WARN)'
"""

if not path.exists('./Log'):
    mkdir('./Log')

class Log:
        """ Define the Logging settings, log file, format, Log on/off
        """
        def __init__(self):
                log_list = ["Log/log_mass.csv", "Log/log_heat_gas.csv","Log/log_any.csv","Log/log_heat_sor.csv"]
                for f in log_list:
                        if path.exists(f):
                                remove(f)

                # Any
                lg_any = logging.getLogger('Any')
                handler_any = logging.FileHandler(filename = "Log/log_any.csv",mode = 'a')
                handler_any.setFormatter(logging.Formatter("%(message)s"))
                lg_any.setLevel(logging.DEBUG) # log on
                #lg_any.setLevel(logging.WARN)   # log off
                lg_any.addHandler(handler_any)

                """set the logger"""
                Log.__instance = self     


        def log_set_any(self,msg):
                log = logging.getLogger('Any')
                log.debug(msg)


""" creating instance of Logger """
l = Log()


""" output any message to "Log/log_any.csv" file """
def log_any_msg(msg):
        l.log_set_any(msg)
        return 0