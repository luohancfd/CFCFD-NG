Transient Recorder system for the T4 shock tunnel
Linux device drivers & debug software
Zane Smith - July 2004

Requires a Linux 2.4.x Kernel

Installation:
1. Make sure FUSD is installed.
2. Run "tranrec" as root.

Running
1. In a new shell, run "python dbdebug.py" as root.

Technical Information:
Upon running, the driver is in the default "Input" mode.
From this mode there are 9 options, initiated by writing a number to /dev/tranrec
2: Arm the trigger units
3: Trigger the trigger units
4: Enter the trigger status read mode
5: Enter the timebase status read mode
6: Enter the A/D card S register read mode
7: Enter the A/D card buffer read mode
8: Enter the direct read mode
9: Enter the direct write mode
11: Enter the Card select mode

Selecting options 4-7 initialise the read mode. Reading from /dev/tranrec will return unprefixed hexadecimal bytes, except in the case of option 7, where it returns the 16kB buffer as ASCII characters.
In order to perform a direct read, first enter option 8, then write the port to be read to /dev/tranrec (in decimal). Reading from /dev/tranrec will then return the hexadecimal byte from that port.
In order to perform a direct write, first enter option 9, then write the port to be written to /dev/tranrec (in decimal), then write the data to be written to /dev/tranrec (in decimal).
In order to select a new card/channel, enter option 11, then write a hexadecimal byte to /dev/tranrec where the high nibble is the channel number, and the low nibble is the card number.

After these actions, the driver returns to input mode.
