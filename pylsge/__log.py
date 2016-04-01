# -*- coding:utf-8
import logging
def initlog(name='defaultlog',filename='default.log'):
    logger = logging.getLogger(name)
    hdlr = logging.FileHandler(filename)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s \n%(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    return logger