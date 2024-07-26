#!/home/ahutko/miniconda3/bin/python

import os
import psycopg2
from datetime import datetime, timedelta

#--------- get picks from PNSN database for various source types
# https://internal.pnsn.org/LOCAL/WikiDocs/index.php/Accessing_the_AQMS_databases

#----- Convert unix/epoc time to truetime which includes leap seconds.
#      Input is a timestamp, e.g. 1568130188.11
#      Output is a datetime object, e.g. datetime.datetime(2019, 9, 10, 15, 42, 41, 117380)
def unix_to_true_time(unixtime):
    leap_seconds = {2272060800: 10 ,2287785600: 11 ,2303683200: 12 ,2335219200: 13 ,2366755200: 14 ,2398291200: 15 ,2429913600: 16 ,2461449600: 17 ,2492985600: 18 ,2524521600: 19 ,2571782400: 20 ,2603318400: 21 ,2634854400: 22 ,2698012800: 23 ,2776982400: 24 ,2840140800: 25 ,2871676800: 26 ,2918937600: 27 ,2950473600: 28 ,2982009600: 29 ,3029443200: 30 ,3076704000: 31 ,3124137600: 32, 3345062400: 33 ,3439756800: 34 ,3550089600: 35 ,3644697600: 36 ,3692217600: 37}
    time1900 = unixtime + 2208988800
    seconds_to_sub = 0
    for utime in leap_seconds:
        if ( time1900 >= utime ):
            seconds_to_sub = leap_seconds[utime] - 10
    t = datetime.utcfromtimestamp(unixtime) - timedelta(seconds=seconds_to_sub)
    return t

#----- Check which evids have already been done
evids_done, evids_tried = [],[]
f = open('event_classification_list.txt')
lines = f.readlines()
for line in lines:
    try:
        if 'Analyst' in line and 'N/A' not in line:
            evids_done.append(int(line.split()[0]))
        else:
            evids_tried.append(int(line.split()[0]))
    except:
        pass

#----- Connect to database

dbname = os.environ['AQMS_DB']
dbuser = os.environ['AQMS_USER']
hostname = os.environ['AQMS_HOST1']  # check which is currently secondary 
dbpass = os.environ['AQMS_PASSWD']

conn = psycopg2.connect(dbname=dbname, user=dbuser, host=hostname, password=dbpass)
cursor = conn.cursor()

now = datetime.now()
two_days_ago = now - timedelta(days=2)
str_two_days_ago = two_days_ago.strftime('%Y-%m-%d %H:%M:%S')

evinfo = {}
cursor.execute('select o.evid, o.orid, o.datetime, o.lat, o.lon, o.depth, o.distance, o.wrms, o.algorithm, o.rflag, n.magnitude, n.magtype, n.uncertainty, n.nsta, n.magalgo from origin o inner join event e on o.evid = e.evid inner join netmag n on n.magid = e.prefmag where ( o.algorithm like (%s) or o.algorithm like (%s) or o.algorithm like (%s) ) and to_timestamp(o.datetime) > (%s) and e.selectflag = (%s) order by o.datetime desc', ( 'BIND%', 'HYP%', 'SUBNE%', str_two_days_ago, 1 ) )
for record in cursor:
    evid = record[0]
    orid = record[1]
    odate = unix_to_true_time(record[2])
    lat = record[3]
    lon = record[4]
    dep = record[5]
    mindist = record[6]
    orms = record[7]
    algorithm = record[8]
    mag = record[10]
    unc = record[11]
    nsta = record[12]
    magalgo = record[13]
    evinfo[evid] = [ orid, odate, lat, lon, dep, mindist, orms, mag, unc, nsta, magalgo, algorithm ]

#----- Now do subnet triggers
cursor.execute('select e.evid, o.orid, o.datetime from event e inner join origin o on e.prefor = o.orid where e.etype = (%s) and to_timestamp(o.datetime) > (%s) and e.subsource = (%s) order by o.datetime desc',( 'st', str_two_days_ago, 'RT1' ) )
for record in cursor:
    evid = record[0]
    orid = record[1]
    odate = unix_to_true_time(record[2])
    lat, lon, dep, mindist, orms, mag, unc, nsta, magalgo, algorithm = 0, 0, 0, 0, 0, 0, 0, 0, 'None', 'None'
    evinfo[evid] = [ orid, odate, lat, lon, dep, mindist, orms, mag, unc, nsta, magalgo, algorithm ]

#----- Now print out evids to run
for evid in evinfo:
    if evid not in evids_done:
        if evid in evids_tried and 'HYP' in evinfo[evid][11]:
            print('')
            print('./run_Automated_Surface_Event_Detection5.py ' + str(evid) )
        elif evid not in evids_tried:
            print('')
            print('./run_Automated_Surface_Event_Detection5.py ' + str(evid) )

