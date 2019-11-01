#! /usr/bin/env python3

# https://www.fullstackpython.com/blog/build-first-slack-bot-python.html
# export SLACK_BOT_TOKEN=xoxb-2984982615-762823711127-RU9n3ojzHVkmqwI0VLVakxMq
# source autoprimer/bin/activate
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
                return message.lower(), event["channel"], event["user"]
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

def handle_command(command, channel, poster):
    """
        Executes bot command if the command is known
    """
    # Default response is help text for the user
    wkr_commands = ['find_primers', 'connection_status', 'hello'] # list of commands needs to be updated as Worker class is
    default_response = "Not sure what you mean. Try *{}*.".format(wkr_commands)
    now = datetime.datetime.now()

    # main loop to look for commands    
    found_command = False
    for com in wkr_commands:
        if command.startswith(com):
            found_command = True

            # parse the command options
            options = parse_command(command, com)

            # debug
            # print(now, command, options, channel, poster)

            # add the job to the schedule
            SCHEDULE[now] = ntp.Worker(now, command, options, channel, poster)
    
    # Default response
    if not found_command:
        catch_message = f"I didn't recognize that command. Try using one of these commands : {','.join(wkr_commands)}"
        post(catch_message, channel)

def post(message, channel):
    slack_client.api_call(
        "chat.postMessage",
        channel=channel,
        text=message
    )

########################
###### Run
########################

# instantiate Slack client
slack_client = SlackClient(os.environ.get('SLACK_BOT_TOKEN'))
# autobot's user ID in Slack: value is assigned after the bot starts up
autobot_id = None

# DEFAULT_FOLDER = '/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autobot_beta/'
DEFAULT_FOLDER = '/mnt/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autobot_beta/'


# constants
RTM_READ_DELAY = 1 # 1 second delay between reading from RTM
# COMMANDS = {
#     "make_primers" : make_primers,
#     "hello": hello,
#     "help" : help,
#     "status" : status,
# }
MENTION_REGEX = "^<@(|[WU].+?)>(.*)"
SCHEDULE = {}

if __name__ == "__main__":
    if slack_client.rtm_connect(with_team_state=False):
        print("AutoBot connected and running!")
        # Read bot's user ID by calling Web API method `auth.test`
        autobot_id = slack_client.api_call("auth.test")["user_id"]
        while True:
            command, channel, poster = parse_bot_commands(slack_client.rtm_read())
            if command:
                handle_command(command, channel, poster)
            time.sleep(RTM_READ_DELAY)

            # post any start messages
            for wokr in SCHEDULE:
                # post any start messages
                if SCHEDULE[wokr].status == 'init':
                    post(SCHEDULE[wokr].start_message(), channel)
                    # run anything in the scheduler
                    post(SCHEDULE[wokr].run(), channel)

                # post any done messages
                if SCHEDULE[wokr].status == 'done':
                    stopmess = post(SCHEDULE[wokr].done_message(), channel)
                    if stopmess: 
                        post(stopmess)
    else:
        print("Connection failed. Exception traceback printed above.")




########################
###### Commands
########################

# def make_primers(options):
#     """
#     Make primers using folder (options)
#     """
#     if not options:
#         options = DEFAULT_FOLDER
#     now = datetime.datetime.now()

#     # add the job to the schedule
#     SCHEDULE[now] = {
#         'bot':Worker(now, command, options),
#         'status' : 'init'
#     }

#     response = f"Starting 'make_primers' at " + now.strftime("%Y-%m-%d.%H-%M")
#     slack_client.api_call(
#         "chat.postMessage",
#         channel=channel,
#         text=response or default_response
#     )
#     # run the command
#     ntp.submit_folder(options)

#     # post another message
#     now = datetime.datetime.now()
#     response = f"Finshed 'make_primers' at " + now.strftime("%Y-%m-%d.%H-%M")
#     slack_client.api_call(
#         "chat.postMessage",
#         channel=channel,
#         text=response or default_response
#             )

# def help(options):
#     now = datetime.datetime.now()
#     response = f"Try using one of these commands : {', '.join(list(COMMANDS.keys()))}"
#     if options:
        


# def status(options):
#     find_folder = os.path.isdir(DEFAULT_FOLDER)
#     if find_folder:
#         response = "Everything looks good, ready to run"
#     else:
#         response = "I can't find the autoprimer gene folder - it is possible that the connection is bad"
#     slack_client.api_call(
#         "chat.postMessage",
#         channel=channel,
#         text=response
#     )    

# def hello(options):
#     response = f"Hello <@{poster}>! Looking forward to designing some primers for you."
#     slack_client.api_call(
#         "chat.postMessage",
#         link_names=1,
#         channel=channel,
#         text=response
#     )    