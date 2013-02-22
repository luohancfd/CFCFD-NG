# fieldsForRunDescription.tcl
# Common definitions for the run description database.
#
# PJ, April 2004
# --------------------------------------------------------------

if { $::debugLevel } {
    puts "sourcing fieldsForRunDescription.tcl"
}

# The following set of (incomplete) regular-expressions identifies the fields,
# assuming that field labels are at the start of their respective lines.

set pattern(facility_name)   [list "Faci"]
set pattern(shot_number)     [list "Run" "Shot" "Num"]
set pattern(basic_file_name) [list "Basic"]
set pattern(date_string)     [list "Date" "Dat"]
set pattern(project)         [list "Proj"]
set pattern(blame)           [list "Bla"]
set pattern(condition)       [list "Cond"]
set pattern(nozzle)          [list "Nozz"]
set pattern(diaphragm1)      [list "Diaphragm\." "Diaphragm " "Diaphragm1" "Dia1"]
set pattern(diaphragm2)      [list "Diaphragm2" "Dia2"]
set pattern(diaphragm3)      [list "Diaphragm3" "Dia3"]
set pattern(reservoir)       [list "Res"]
set pattern(driver0)         [list "Driv"]
set pattern(driver1)         [list "1st_Dr"]
set pattern(driver2)         [list "2nd_Dr"]
set pattern(shocktube)       [list "Shoc"]
set pattern(acceltube)       [list "Accel"]
set pattern(dumptank)        [list "Evac"]
set pattern(jackoff)         [list "Jac"]
set pattern(piston)          [list "Pist"]
set pattern(notes)           [list "AllNotes"]

# A matching list of field names for the database
set fieldNames [list facility_name shot_number basic_file_name date_string \
		    project blame condition nozzle \
		    diaphragm1 diaphragm2 diaphragm3 reservoir \
		    driver0 driver1 driver2 \
		    shocktube acceltube dumptank \
		    jackoff piston notes]

set nullString {\N}

set sqlType(facility_name)   "text"
set sqlType(shot_number)     "integer"
set sqlType(basic_file_name) "text"
set sqlType(date_string)     "text"
set sqlType(project)         "text"
set sqlType(blame)           "text"
set sqlType(condition)       "text"
set sqlType(nozzle)          "text"
set sqlType(diaphragm1)      "text"
set sqlType(diaphragm2)      "text"
set sqlType(diaphragm3)      "text"
set sqlType(reservoir)       "text"
set sqlType(driver0)         "text"
set sqlType(driver1)         "text"
set sqlType(driver2)         "text"
set sqlType(shocktube)       "text"
set sqlType(acceltube)       "text"
set sqlType(dumptank)        "text"
set sqlType(jackoff)         "text"
set sqlType(piston)          "text"
set sqlType(notes)           "text"

# ----------------------------------------------------------------------------
