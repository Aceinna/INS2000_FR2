#### Python script to log data

Prerequisite: Python 3.x

- Femtomes INS2000_FRII data logging: femtomes_logging.py
- SJ Novatel GNSS receiver + iMAR IMU data logging:  novatel_imar_logging.py

All you need to do is to change the serial port number `ser = serial.Serial('com140',230400,parity='N',bytesize=8,stopbits=1,timeout=None) `

in this line, and start logging



#### Decoder executable to parse data

Change the span_decoder.exe1 to span_decoder.exe, and run in the cmd line

```
span_decoder.exe <data.bin> 1
```

it works for SJ NovAtel system data and Femtomes data, and decodes the bin files to csv files. Among the decoded csv files, "_gnss.txt" and "_ins.txt" are GNSS and INS solutions respectively. 

