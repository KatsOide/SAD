if {![package vsatisfies [package provide Tcl] 8.5]} { return }
if {[package vcompare [package provide Tcl] 8.5a6] != 0} { return }
package ifneeded Tk 8.5a6 [list load [file join $dir .. libtk8.5.dylib] Tk]
