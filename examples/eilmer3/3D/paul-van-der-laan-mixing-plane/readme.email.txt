From: Paul van der Laan [mailto:paulvdlaan@gmail.com]
Sent: Tue 12/1/2009 3:50 PM
To: Paul Petrie-Repar
Subject: Mixing plane files
 
Paul,

The mixing plane files as you requested.

Keep in mind that a simulation without omega_rotor = 2000 needs the
following changes:

- downstream mp lua files: omega_z = 0 instead of 2000
- use testcase_full_run.sh (not testcase_full_run_omega.sh)
- omega_z = 0.0 instead of 2000.0 in testcase.py

Greetings,
Paul 

