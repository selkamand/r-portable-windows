# -*- tcl -*-
# Tcl package index file, version 1.1
#
# Note sqlite*3* init specifically
#
if {[package vsatisfies [package provide Tcl] 9.0-]} {
    package ifneeded sqlite3 3.40.0 \
	    [list load [file join $dir tcl9sqlite3400.dll] Sqlite3]
} else {
    package ifneeded sqlite3 3.40.0 \
	    [list load [file join $dir sqlite3400.dll] Sqlite3]
}
