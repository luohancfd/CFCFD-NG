proc wait_a_moment {} {
    # Wait for a little while still processing events.
    after 200 [list set done 1]
    vwait done
}; # end proc
