#!/usr/bin/env python
"""
e3_load_balance.py -- Program to do some static load balancing.

The idea of this program routines is to re-order a supplied
hostfile in such a way as to give a decent load balancing when
running MPI in an oversubscribed manner. What we are trying to solve
is what is referred to as the "multiprocessor scheduling
problem", which may be stated as [paraphrasing Graham (1969)]:

Let us suppose we are given 'n' identical processing
units P_i, i = 1, ..., n, and a set of tasks,
T = {T_i, ...., T_r} which is to be processed by
the P_i. We also have a measure of the time it
takes to complete each task, mu(T_i).  The total
time taken to complete the tasks is 'w'.

In our restricted problem, we wish to find
the best way to assign the T_i to the P_i
which minimises 'w'. We have no priority
list for the task, that is, any task may 
start at any time the resources are free,
and is not conditional on a prior task
finishing.

We will do the task assignment in a static manner. The
estimates of computed load per block are simply based
on the block dimensions. These estimates may be poor
if we are doing a lot of work in certain cells compared
to others due to chemistry or thermal nonequilibrium,
for example.

The algorithm chosen [Graham (1969)] is not provably
optimal, but it should generally give a good result
for a relatively cheap computation. It is possible to
find the optimimum task packing by brute force enumerating
all possible ways the tasks could be packed into the
processing units, but the computational expense of
this increases exponentially with number of tasks
and number of processing units.

Reference:
---------
Graham, R. L. (1979)
Bounds on Multiprocessing Timing Anomalies,
SIAM Journal of Applied Mathematics, 17:2, pp. 416--429

.. Author: Rowan J. Gollan
"""

import sys
import os
import ConfigParser
from copy import copy
from optparse import OptionParser

OUT_EXT = "mpimap"

parser = OptionParser()
parser.add_option("-j", "--job", dest="jobname", type="string",
                  help="base file name of job",
                  metavar="JOBNAME")
parser.add_option("-i", "--input-hostfile", dest="hosts_in", type="string",
                  help="supplied hostfile (input)",
                  metavar="HOSTS_IN")
parser.add_option("-n", "--number-of-procs", dest="nprocs", type="int",
                  help="number of processors to use",
                  metavar="NP")
parser.add_option("-s", "--sweep-range", dest="sweep", type="string",
                  help="sweep through range of 'lower_np:upper_np' testing load balance quality")
#parser.add_option("-o", "--output-hostfile", dest="hosts_out", type="string",
#                  help="output hostfile with entries re-ordered",
#                  metavar="HOSTS_OUT")

def create_task_list(cfg_fname):
    """
    Return a list of tuples of form (task_id, task_load).

    The returned list will have one entry per block.
    It will also be sorted based on load grom highest
    to lowest.
    
    :param ctrl_fname: Name of eilmer control file.
    """
    cp = ConfigParser.ConfigParser()
    cp.read(cfg_fname)
    nblks = int(cp.get("global_data", "nblock"))

    tasks = []
    total_load = 0
    for ib in range(nblks):
        sec_name = "block/" + str(ib)
        nni = int(cp.get(sec_name, "nni"))
        nnj = int(cp.get(sec_name, "nnj"))
        nnk = int(cp.get(sec_name, "nnk"))
        # Estimate task load as: nni*nnj*nnk
        load = nni*nnj*nnk
        tasks.append((ib, load))
        total_load += load
    #
    # Sort tasks by load (2nd element in tuple)
    tasks.sort(key=lambda x: x[1], reverse=True)
    return tasks, total_load

def load_balance(tasks, nprocs):
    procs = []
    for i in range(nprocs): procs.append([])
    proc_loads = [0,]*nprocs
    task_hosts = []
    for (task_id, task_load) in tasks:
        # The algorithm is simply this:
        # Given a list sorted from heaviest load to lightest,
        # at any stage, assign task to the processing unit
        # with the minimum load currently.
        min_load = min(proc_loads)
        imin = proc_loads.index(min_load)
        proc_loads[imin] += task_load
        procs[imin].append(task_id)
        task_hosts.append(imin)
    return procs, proc_loads, task_hosts

def main():
    (options, args) = parser.parse_args()
    
    if options.jobname is None:
        print "The base file name for the job must be specified with"
        print "--job=JOBNAME"
        print parser.print_help()
        sys.exit(1)


    cfg_fname = options.jobname + ".config"
    tasks, total_load = create_task_list(cfg_fname)

    if options.nprocs:
        nprocs = options.nprocs
    elif options.hosts_in:
        hosts = parse_hosts_file(options.hosts_in)
        nprocs = len(hosts)
    else:
        print "The number of processors has not been supplied."
        print "This should be supplied directly with: -np NP"
        print "or indirectly by giving a hostsfile to be parsed: --input-hostfile=HOSTS"
        print "Bailing out!"
        print parser.print_help()
        sys.exit(1)

    if nprocs > len(tasks):
        print "e3loadbalance.py: Quitting without doing anything"
        print "because nprocs > nblocks."
        print "nprocs= ", nprocs
        print "nblocks= ", len(tasks)
        print "Done."
        sys.exit(0)
    #
    
    procs, proc_loads, task_hosts = load_balance(tasks, nprocs)
    #
    # Now write output file in INI format.
    fname = options.jobname + "." + OUT_EXT
    print "Writing out mpi rank-to-block map: ", fname
    f = open(fname, 'w')
    f.write("[global]\n")
    f.write("nrank = %d\n" % nprocs)
    for ip, tsk_list in enumerate(procs):
        f.write("[rank/%d]\n" % ip)
        f.write("nblock = %d\n" % len(tsk_list))
        f.write("blocks = ")
        for t in sorted(tsk_list):
            f.write("%d " % t)
        f.write("\n")
    #
    f.close()
    print "Done."
    #
    if options.sweep:
        tks = options.sweep.split(":")
        lower_n = int(tks[0])
        upper_n = int(tks[1])
        if upper_n < lower_n:
            print "Error in sweep range string. upper_n should be greater than lower_n."
            print "lower_n= ", lower_n
            print "upper_n= ", upper_n
            print "Bailing out!"
            sys.exit(1)
        print "Performing sweep of nprocs."
        f = open('load-balance.dat', 'w')
        f.write("# nprocs    packing-quality  speedup\n")
        for n in range(lower_n, upper_n+1):
            if n > len(tasks):
                break
            print "nprocs= ", n
            procs, proc_loads, task_hosts = load_balance(tasks, n)
            max_load = max(proc_loads)
            min_load = min(proc_loads)
            delta_cells = max_load - min_load
            pq = 1.0 - float(max_load - min_load) / max_load
            speedup = float(total_load) / max_load
            f.write("%03d  %d %12.6f  %12.6f\n" % (n, delta_cells, pq, speedup))
        f.close()
        print "Written data to: load-balance.dat"
    
    print "e3loadbalance.py: Done."
    return

if __name__ == '__main__':
    main()


    





