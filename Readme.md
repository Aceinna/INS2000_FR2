### Python script to log data

#### Prerequisite: 

​	Python 3.x

#### scripts:

- Femtomes INS2000_FRII data logging: femtomes_logging.py
- San Jose (SJ) system:  Novatel GNSS receiver + iMAR IMU data logging:  novatel_imar_logging.py

  1.  change the serial port number in the line 

     ```python
     ser = serial.Serial('com140',230400,parity='N',bytesize=8,stopbits=1,timeout=None) 
     ```

  2. change the lever arm parameters (refer to respective product manual)

     - Femtomes system: 

       ```shell
       SETINSTRANSLATION ANT1 <parameters>
       SETINSTRANSLATION DUALANT <parameters>
       SETINSROTATION RBV <parameters> # not necessary if aligned according to manual
       ```

     - NovAtel system: 

       ```shell
       setimutoantoffset <parameters>
       ```

       

### Decoder executable to parse data

Change the span_decoder.exe1 to span_decoder.exe, and run in the cmd line

```
aceinna_decoder.exe <data.bin> 1
```

it works for SJ NovAtel system data and Femtomes data, and decodes the bin files to csv files. Among the decoded csv files, "_gnss.txt" and "_ins.txt" are GNSS and INS solutions respectively. 

