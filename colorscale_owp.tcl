proc colorscale_owp { {reverse 0} {count 0} } {
  display update off
  set mincolorid [expr [colorinfo num] - 1]
  set maxcolorid [expr [colorinfo max] - 1]
  set colrange [expr $maxcolorid - $mincolorid]
  set colhalf [expr $colrange / 2]

  for {set i $mincolorid} {$i <= $maxcolorid} {incr i} {
    set x [expr ($i - $mincolorid) / double($colrange)]

    if { $count != 0 } {
      set nx [expr {int($count * $x) / double($count)}]
      set x $nx
    }

    if {$x > 1.0} {
      set x 1.0;
    }

    if {$x < 0.0} {
      set x 0.0;
    } 

    set r 0.0;
    if {$x <= 0.5} {
      set r 1.0;
    } elseif {$x > 0.5} {
      set r [expr 1.5 - $x];
    }
    
    set g 0.0;
    if {$x <= 0.5} {
      set g [expr 0.5 + $x];
    } elseif {$x > 0.5} {
      set g [expr 2.0 - 2 * $x];
    }

    set b 0.0;
    if {$x < 0.5} {
      set b [expr 2.0 * $x];
    } elseif {$x >= 0.5} {
      set b 1.0;
    }

    if { $reverse } {
      color change rgb [expr $mincolorid + ($maxcolorid - $i)] $r $g $b
    } else {
      color change rgb $i $r $g $b
    }
  }

  display update ui
  display update on
}
