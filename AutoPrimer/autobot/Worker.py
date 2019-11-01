#! /usr/bin/env python3

import datetime
import AutoPrimer as ntp
import os

class Worker(object):

    def __init__(self, name, command, options, channel, poster):
        # date time stamp from scheduler
        self.name = name
        self.status = 'init' # init, running, done, expired
        
        # which command this worker should execute
        self.command = command
        
        # what options to pipe to the command
        self.options = options

        self.channel = channel
        self.poster = poster

        # defines possible commands
        self.comms = {
            'find_primers' : self.find_primers,
            'connection_status' : self.connection_status,
            'hello' : self.hello
        }

        # all commands from comms should go into one of these two categories
        self.message_commands = ['connection_status', 'hello']
        self.action_commands = ['find_primers']


    def start_message(self):
        # start and stop message commands
        if self.command in self.action_commands:
            now = datetime.datetime.now()
            mess = f"Starting {self.command} at " + now.strftime("%Y-%m-%d.%H-%M")
            return mess
        
        # start only message commands
        if self.command in self.message_commands:
            return self.comms[self.command]()

    def run(self):
        if self.command in self.action_commands:
            self.comms[self.command]()

    def done_message(self):
        if self.command in self.action_commands:
            now = datetime.datetime.now()
            mess = f"Finshed {self.command} main at " + now.strftime("%Y-%m-%d.%H-%M")
            common_dict[self.name]['status'] = 'closed'
            self.status = 'expired'
            return mess
    
    #######
    ##### Commands
    #######

    def find_primers(self):
        """
        Initiates AutoPrimer
        """
        self.status = 'running'
        if not self.options:
            self.options = '/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer'
        find_folder = os.path.isdir('/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer')
        if find_folder:
            ntp.submit_folder(self.options)
            self.status = 'done'
        else:
            self.status = 'failed'
    
    def connection_status(self):
        """
        Checks the connection status to a given folder
        """
        if self.options:
            folder = self.options
            find_folder = os.path.isdir(self.options)
            
        else:
            folder = '/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer'
            find_folder = os.path.isdir('/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer')
        
        if find_folder:
            response = f"Connection established with {folder}. Everything looks good, ready to run"
        else:
            response = f"I can't find the {folder} - it is possible that the connection is bad"
        
        self.status = 'expired'
        return response

    def hello(self):
        """
        Replies to the user with a friendly hello message.
        """
        response = f"Hello <@{self.poster}>! Looking forward to designing some primers for you."
        self.status = 'expired'
        return response