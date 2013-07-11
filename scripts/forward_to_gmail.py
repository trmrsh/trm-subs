#!/usr/bin/env python

"""
Script to forward email from the Warwick IMAP server's INBOX to googlemail.
The script stores the timestamps of the forwarded messages (to the
nearest second) and uses these to avoid multiple forwarding of the
same message. Probably some better way but it should work OK, with the
slight risk that some messages are missed because they arrived at the
same time as another already-forwarded one.

Arguments:

wuser   -- Warwick account user name
wpass   -- Warwick account password
guser   -- googlemail account user name
gpass   -- googlemail account password

"""
#

import os
import sys
import mimetypes
import email.Parser
import email.Utils
import imaplib
import smtplib
import exceptions
import time

def search_date(timestamp):
    """Converts time stamp into format needed for IMAP search"""
    dtup = time.localtime(timestamp)
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    return str(dtup[2]) + '-' + months[dtup[1]-1] + '-' + str(dtup[0])

def error_handler(message):
    imap.close()
    imap.logout()
    raise SubsError(message)

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print 'usage: wuser wpass guser gpass'
        exit(0)

    (wuser,wpass,guser,gpass) = sys.argv[1:]

    home = os.getenv('HOME')
    if home is None:
        raise SubsError('Could not find home directory')

    # Maximum number of timestamps to save
    MAXSAVE = 1000

    # try to read timestamps of most recently forwarded messages
    tfile = os.path.join(home, '.forward_to_gmail')
    try:
        fptr = open(tfile)
        most_recent = []
        for line in fptr:
            most_recent.append(int(line))
        fptr.close()
    except IOError:
        most_recent = []
        
    # login to the Warwick imap server
    imap = imaplib.IMAP4_SSL('myimapmail.warwick.ac.uk', 993)
    imap.login(wuser, wpass)

    # select the inbox
    imap.select()

    # search for messages, starting from 2 days ago
    stime = int(time.time()) - 2*86400

    start_date = search_date(stime)
    print 'Searching for mail sent after ' + start_date
    ret, data = imap.search(None, 'SINCE', start_date)
    if ret != 'OK': error_handler('imap.search return = ' + ret)

    last_time = most_recent
    if data[0] != '':

        # we have found some messages
        print 'Found ' + str(len(data[0].split())) + ' messages.'

        # create email parser
        parser = email.Parser.Parser()

        # connect to the gmail smtp server
        smtps = smtplib.SMTP('smtp.gmail.com', 587)
#        smtps.set_debuglevel(1)
        smtps.ehlo()
        smtps.starttls()
        smtps.ehlo()
        smtps.login(guser, gpass)

        nfor = 0 
        for num in data[0].split():

            # BODY.PEEK as opposed to BODY stops the mail being marked as read
            # Just get the date to start with to reduce the amount of download
            # since the majority will have already been sent.
            ret, data = imap.fetch(num, '(BODY.PEEK[HEADER.FIELDS (DATE)])')
            if ret != 'OK': 
                smtps.close()
                error_handler('imap.fetch of date returned = ' + ret)

            date      = data[0][1][6:].strip()
            datetime  = email.Utils.parsedate(date)
            if datetime is not None:
                ftime = int(time.mktime(datetime))

                # Only bother with those we have not already marked as done. This is
                # done by time (to the nearest second)
                if ftime not in most_recent:
                    
                    print 'Fetching message',num,'...'
                    ret, data = imap.fetch(num, '(BODY.PEEK[])')
                    if ret != 'OK': 
                        smtps.close()
                        error_handler('imap.fetch returned = ' + ret)

                    msg = parser.parsestr(data[0][1])

                    # send off to gmail, avoiding those marked as spam
                    if not msg['Subject'].startswith('***SPAM***'):
                        try:
                            smtps.sendmail(msg['From'], 'tom.r.marsh@gmail.com', msg.as_string())
                            print 'Forwarded message',num,'from',msg['From'],'to googlemail'
                            nfor += 1
                            # wait a bit before sending next one off to avoid overloading the server
                            time.sleep(2)

                        except smtplib.SMTPDataError, err:
                            print 'Failed to forward message',num,'from',msg['From'],'with error: ' + str(err)
                    
                    # even if we fail, we store the time stamp to avoid repeated attempts
                    # to forward hard cases.
                    most_recent.append(ftime)
            
            else:
                print 'Datetime tuple from message ' + str(num) + ' was None; this may fix itself later.'

        smtps.close()

        # store MAXSAVE most recent timestamps. So long as I don't get this many every 2 days,
        # things should be OK ....
        most_recent.sort()
        fptr = open(tfile, 'w')
        for tstamp in most_recent[-MAXSAVE:]:
            fptr.write(str(tstamp) + '\n')
        fptr.close()

    
    imap.close()
    imap.logout()

    if nfor == 1:
        print '1 message was forwarded.'
    else:
        print nfor,'messages were forwarded.'

    exit(0)
        


