#### Python script to logging data

Prerequisite: Python 3.x

- Femtomes INS2000_FRII data logging: femtomes_logging.py
- SJ Novatel GNSS receiver + iMAR IMU data logging:  novatel_imar_logging.py

All you need to do is to change the serial port number `ser = serial.Serial('com140',230400,parity='N',bytesize=8,stopbits=1,timeout=None) `

in this line, and start logging