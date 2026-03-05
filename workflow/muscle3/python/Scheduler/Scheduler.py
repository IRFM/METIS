import logging
import os
import math

import numpy as np
import matplotlib.pyplot as plt


from libmuscle import Grid, Instance, Message
from ymmsl import Operator

def Scheduler() -> None:
    """Just to provide time evolution to METIS"""

    logger = logging.getLogger()
    logger.info('time manager started')
    
   # instance = Instance({
   #         Operator.O_F: ['next_time'],   
   #         Operator.F_INIT: ['status']})

    instance = Instance({
            Operator.S: ['next_time'],   
            Operator.O_I: ['status']})
            
    logger.info('Read parameters')
  
    t_start = instance.get_setting('time_start', 'float')
    t_stop = instance.get_setting('time_stop', 'float')
    dt = instance.get_setting('time_step', 'float')
    
    logger.info('Start @ %f s',t_start)
    logger.info('Stop @ %f s',t_stop)
    logger.info('Time step = %f s',dt)
    
    t_cur  = t_start
    t_next = t_start
        
    loop=0
    total_loop=0
      
    while instance.reuse_instance():
 
        logger.info('In reuse loop')

        while t_next < t_stop:
     

                
            logger.info('--------------------------')
            logger.info('current time = %f s',t_cur)   
                            

            # next time step
            t_next = t_cur + dt

            logger.info('step time for the next loop: %f s, next time to compute sent : %f s',dt,t_next)
            
            out_msg = Message(t_cur, t_next, t_cur)
            instance.send('next_time', out_msg)
            logger.info('message sent')


             # not yet used for initialisation with steady time
            t_cur = t_next

            # synchronisation with METIS
            in_msg = instance.receive('status')
            logger.info('message received')
           	
if __name__ == '__main__':
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    Scheduler()
