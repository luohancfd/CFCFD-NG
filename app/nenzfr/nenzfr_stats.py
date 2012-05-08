"""
nenzfr_stats.py -- Flow statistics function needed by the main program.
"""

import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------

def print_stats(sliceFileName,jobName,coreRfraction):
    """
    Display statistics of flow properties at the nozzle exit.
    """
    print "Nozzle-exit statistics:"
    fp = open(sliceFileName, 'r')
    # Keep a list of variables in order of appearance.
    varLine = fp.readline().strip()
    items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # print "variable_list=", variable_list
    # Store the data in lists against these names.
    data = {}
    for var in variable_list:
        data[var] = []
    for line in fp.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        assert len(items) == len(variable_list)
        for i in range(len(items)):
            data[variable_list[i]].append(float(items[i]))
    fp.close()
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    y_edge = ys[-1] * coreRfraction
    #
    # Compute and print area-weighted-average core flow values.
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.stats','w')
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
    for var in variable_list:
        if var in exclude_list: continue
        A = 0.0; F = 0.0;
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            if j == 0:
                y0 = 0.0
            else:
                y0 = 0.5*(ys[j-1]+ys[j])
            y1 = 0.5*(ys[j]+ys[j+1])
            dA = y1**2 - y0**2
            F += data[var][j] * dA
            A += dA
        mean = F/A
        # Identify low and high values.
        diff_minus = 0.0
        diff_plus = 0.0
        count = 0.0
        stddev = 0.0
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            diff = data[var][j] - mean
            diff_minus = min(diff, diff_minus)
            diff_plus = max(diff, diff_plus)
            count += 1
            stddev += diff**2
        # Calculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%10s  %12.4g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.4g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #    
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return
