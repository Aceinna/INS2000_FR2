#!/usr/bin/python
import serial
import math
import time
import datetime

def getweeknum(weekseconds):
    return math.floor(weekseconds/(7*24*3600))

def getweeksec(weekseconds):
    return weekseconds - getweeknum(weekseconds)*(7*24*3600)

def yearfour(year):
    if year<=80:
        year += 2000
    elif year<1990 and year>80:
        year += 1900
    return year

def isleapyear(year):
    return (yearfour(year)%4==0 and yearfour(year)%100!=0) or yearfour(year)%400==0
                 
def timefromGPS(weeknum,weeksec):
    year = 0
    month = 0
    day = 0
    hour = 0
    minute = 0
    second = 0
    doy = 0
    daypermon = [31,28,31,30,31,30,31,31,30,31,30,31]

    weeknum += getweeknum(weeksec)
    weeksec  = getweeksec(weeksec)
    
    weekmin  = math.floor(weeksec/60.0)
    second   = weeksec - weekmin*60.0
    weekhour = math.floor(weekmin/60)
    minute   = weekmin - weekhour*60
    weekday  = math.floor(weekhour/24)
    hour     = weekhour - weekday*24

    totalday = weekday+weeknum*7
    if totalday<360:
        year = 1980
    else:
        year = 1981
        totalday -= 360
        while True:
            if totalday<365:
                break
            if isleapyear(year): totalday -= 1
            totalday -= 365
            year += 1
    doy = totalday

    if totalday <= daypermon[0]:
        month = 1
    else:
        totalday -= daypermon[0];
        if isleapyear(year): totalday -= 1
        month = 2
        while True:
            if totalday<=daypermon[month-1]:
                break
            else:
                totalday -= daypermon[month-1]
                month += 1
    if month==2 and isleapyear(year): totalday += 1
    day = totalday
    return [year,month,day,hour,minute,second,doy]


def configFemtomes(ser):

    setupcommands7  = ['unlogall\r',\
                 # change the ntrip server config below
                'NTRIPCONFIG NCOM1 client V1 vrs.sixents.com:8002 RTCM32_GNSS_RTK ke010 etav5s\r',\
                'SETINSPROFILE LAND_BASIC\r',\
                # 'alignmentmode automatic\r',\
                # change the lever arm and rotation angle below, according to your installation
                'SETINSTRANSLATION ANT1 -0.28 1.43 0.88 0.20 0.20 0.20\r',\
                'SETINSTRANSLATION DUALANT 0.7 0.0 0.0 0.20 0.20 0.20\r',\
                'SETINSROTATION RBV 0.0 0.0 0.0 0.5 0.5 1.0\r',\
                # 'SETINSTRANSLATION USER 0.1 1.0 0.0 0.20 0.20 0.20\r',\
                'INSCOMMAND ENABLE\r',\
                'LOG RANGECMPB ONTIME 1\r',\
                'LOG RAWEPHEMB ONCHANGED\r',\
                'LOG GLOEPHEMERISB ONCHANGED\r',\
                # 'log GALFNAVEPHEMERISB ONCHANGED\r',\
                # 'log GALINAVEPHEMERISB ONCHANGED\r',\
                'LOG GALEPHEMERIS ONCHANGED\r',\
                'LOG BDSEPHEMERISB ONCHANGED\r',\
                'LOG QZSSEPHEMERISB ONCHANGED\r',\
                'LOG INSCONFIGB ONCHANGED\r',\
                'LOG RAWIMUSXB ONNEW\r',\
                'LOG versionb once\r',\
                'LOG rxstatusb once\r',\
                'LOG inspvaxb ontime 0.1\r',\
                #'log bestposb ontime 0.1\r',\
                'LOG bestgnssposb ontime 0.1\r',\
                'LOG bestgnssvelb ontime 0.1\r',\
                'LOG heading2 onchanged\r',\
	'LOG gpgga ontime 0.2\r',\
                'LOG heading ontime 0.2\r',\
                'LOG NCOM1 gpgga ontime 1\r',\
                'SAVECONFIG\r']

    for cmd in setupcommands7:
        ser.write(cmd.encode())    
 

ser = serial.Serial('com34',460800,parity='N',bytesize=8,stopbits=1,timeout=None)
fname = './data/femtomes_ins330-'
ser.flushInput()

fmode = 'wb'

while True:
    if ser.isOpen(): break

print ('\Port is open now\n')
configFemtomes(ser)
ser.flushInput()

fname += time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()) + '.bin'

with open(fname,fmode) as outf:
    while True:
        try:
            line = ser.readline()
            outf.write(bytes(line))  #line.decode('utf-8')
            
        except:
            #break
            pass

    outf.close()
