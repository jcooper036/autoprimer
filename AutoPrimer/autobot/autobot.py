#! /usr/bin/env python3

#https://www.fullstackpython.com/blog/build-first-slack-bot-python.html

import AutoPrimer as ntp
import os
import time
import datetime
import re
from slackclient import SlackClient

def parse_bot_commands(slack_events):
    """
        Parses a list of events coming from the Slack RTM API to find bot commands.
        If a bot command is found, this function returns a tuple of command and channel.
        If its not found, then this function returns None, None.
    """
    for event in slack_events:
        if event["type"] == "message" and not "subtype" in event:
            user_id, message = parse_direct_mention(event["text"])
            if user_id == autobot_id:
                return message, event["channel"], event["user"]
    return None, None, None

def parse_direct_mention(message_text):
    """
        Finds a direct mention (a mention that is at the beginning) in message text
        and returns the user ID which was mentioned. If there is no direct mention, returns None
    """
    matches = re.search(MENTION_REGEX, message_text)
    # the first group contains the username, the second group contains the remaining message
    return (matches.group(1), matches.group(2).strip()) if matches else (None, None)

def parse_command(command, com):
    """
        Splits the command into 
    """
    options = command.split(com)[1]
    options = options.strip()
    if not options:
        return False
    return options

def handle_command(command, channel):
    """
        Executes bot command if the command is known
    """
    # Default response is help text for the user
    default_response = "Not sure what you mean. Try *{}*.".format(COMMANDS)

    # main loop to look for commands    
    found_command = False
    for com in COMMANDS:
        if command.startswith(com):
            found_command = True
            # parse the command options
            options = parse_command(command, com)
            # run the command
            COMMANDS[com](options)
    
    # Default response
    if not found_command:
        options = True
        COMMANDS['help'](options)

########################
###### Commands
########################

def make_primers(options):
    """
    Make primers using folder (options)
    """
    if not options:
        options = '/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer/'
    now = datetime.datetime.now()
    response = f"Starting '{com}' at " + now.strftime("%Y-%m-%d.%H-%M")
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=response or default_response
    )
    # run the command
    ntp.submit_folder(options)

    # post another message
    now = datetime.datetime.now()
    response = f"Finshed '{com}' at " + now.strftime("%Y-%m-%d.%H-%M")
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=response or default_response
            )

def help(options):
    now = datetime.datetime.now()
    response = f"Try using one of these commands : {', '.join(list(COMMANDS.keys()))}"
    if options:
        response = "I didn't recognize that command. " + response
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=response or default_response
    )

def status(options):
    find_folder = os.path.isdir('/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autobot_beta/')
    if find_folder:
        response = "Everything looks good, ready to run"
    else:
        response = "I can't find the autoprimer gene folder - it is possible that the connection is bad"
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=response
    )    

def hello(options):
    response = f"Hello <@{poster}>! Looking forward to designing some primers for you."
    slack_client.api_call(
        "chat.postMessage",
        link_names=1,
        channel=channel,
        text=response
    )    

########################
###### Run
########################

# instantiate Slack client
slack_client = SlackClient(os.environ.get('SLACK_BOT_TOKEN'))
# autobot's user ID in Slack: value is assigned after the bot starts up
autobot_id = None

# constants
RTM_READ_DELAY = 1 # 1 second delay between reading from RTM
COMMANDS = {
    "make_primers" : make_primers,
    "hello": hello,
    "help" : help,
    "status" : status,
}
MENTION_REGEX = "^<@(|[WU].+?)>(.*)"


if __name__ == "__main__":
    if slack_client.rtm_connect(with_team_state=False):
        print("AutoBot connected and running!")
        # Read bot's user ID by calling Web API method `auth.test`
        autobot_id = slack_client.api_call("auth.test")["user_id"]
        while True:
            command, channel, poster = parse_bot_commands(slack_client.rtm_read())
            if command:
                handle_command(command, channel)
            time.sleep(RTM_READ_DELAY)
    else:
        print("Connection failed. Exception traceback printed above.")


