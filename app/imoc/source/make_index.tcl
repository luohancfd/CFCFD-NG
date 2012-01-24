# make_index.tcl
# Generate the index file for the IMOC package.
# When installed, the shared-object and the tcl scripts will
# all reside in the one directory.

set suffix [info sharedlibextension]

pkg_mkIndex [pwd] \
    imoc$suffix moc_kernel.tcl \
    moc_gui.tcl moc_placard.tcl moc_menu.tcl moc_plot.tcl \
    moc_nodelist.tcl moc_scales.tcl moc_unitproc.tcl \
    moc_syn_cmds.tcl
