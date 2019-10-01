#! /usr/bin/env python3

import datetime
import AutoPrimer as ntp

class Worker(object):

    def __init__(self, name, command, options, channel, poster):
        # date time stamp from scheduler
        self.name = name
        
        # which command this worker should execute
        self.command = command
        
        # what options to pipe to the command
        self.options = options

        self.channel = channel
        self.poster = poster

        # defines possible commands
        self.comms = {
            'find_primers' : self.find_primers()
        }
    
    def start_messge(self):
        # start and stop message commands
        if self.command in ['find_primers']:
            now = datetime.datetime.now()
            mess = f"Starting {self.command} at " + now.strftime("%Y-%m-%d.%H-%M")
            return mess
        # start only message commands
        if self.command in ['status', 'hello']:
            return self.comms[self.command]()

    def run(self, common_dict):
        if self.command in ['find_primers']:
            common_dict[self.name]['status'] = 'running'
            self.comms[self.command]()
            common_dict[self.name]['status'] = 'done'

    def done_messge(self):
        if self.command in ['find_primers']:
            now = datetime.datetime.now()
            mess = f"Finshed {self.command} main at " + now.strftime("%Y-%m-%d.%H-%M")
            common_dict[self.name]['status'] = 'closed'
            return mess
        else:
            return False
    
    #######
    ##### Commands
    #######

    def find_primers(self):
        ntp.submit_folder(self.options)
        return

    def status(self):
        find_folder = os.path.isdir(self.options)
        if find_folder:
            response = "Everything looks good, ready to run"
        else:
            response = "I can't find the autoprimer gene folder - it is possible that the connection is bad"
        slack_client.api_call(
            "chat.postMessage",
            channel=channel,
            text=response
    )    
    
    def hello(self):
        response = f"Hello <@{self.poster}>! Looking forward to designing some primers for you."
        return response