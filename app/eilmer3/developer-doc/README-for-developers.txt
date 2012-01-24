Eilmer3 Developer Notes
=======================
:Author: Peter_Jacobs
:Email: peterj@mech.uq.edu.au
:Date: 04-October-2009


The bzr repository
------------------
.Some suggestions:
- See Rowan's very good notes on using the bzr revision control system.
  These are found in the file 
  ~/cfcfd2/doc/developer-notes/local-development-with-bazaar.html
- The main repository for code developers is the bzr
  repository on git@triton:/archive1/git/cfcfd2/
- On PJ's workstation the parent branch is thus
  bzr+ssh://git@triton/archive1/git/cfcfd2/
- Keep the history of the main repository linear.
  This can be achieved in a couple of ways:
  . Do a pull from the main repository before any commit then,
    so long as no one else has gotten in before you, do a push.
    Your latest commit will then the be the top of the pile.
  . If you have been doing local commits and have accumulated some
    history of revisions that are not in the main repository, 
    rebase your history on that of the main repository.
    Again, your revisions should appear at the top of the pile so
    you can push them, keeping the history in the main repository
    linear. 
- Again, see Rowan's notes on using bzr.
  There is a much more complete discussion of branches and rebase
  in those notes. 
- Don't push broken code (without warning everyone).
  The code should pass the automated test cases discussed below.
  These functional tests will not guarantee working code but they should
  indicate gross errors that are likely to affect many of the code users.


The build
---------
.Tips for minimizing pain:
- Eliminate all compiler warnings.
  It's difficult to spot significant problems in a lot of noise.
- Every so often use Valgrind to check that we have no obvious 
  memory problems.


The automated test cases
------------------------
We're going to need these as more people are fiddling with the code.
In a few example subdirectories, there is a .test script which uses the tcltest
module to automate running of the programs and looking at some of
the output data to see if it has turned out as expected.

.The current "smoke test" is to:
- compile from a clean tree, 
- rsync the examples/eilmer3/ directory and its contents to a convenient
work directory, and then
- run the top-level test script (which simply runs the individual test scripts
located in the various example subdirectories).

[source,shell]
cd ~/cfcfd2/app/eilmer3/build; make TARGET=for_openmpi install
rsync -av ~/cfcfd2/examples/eilmer3 ~/work/
cd ~/work/eilmer3/; ./eilmer3-test.tcl

The cases will take some minutes each but you should see progress
as each test is run.
If things are looking promising, take a break and go get a cup of tea or coffee.
After a few minutes and if all is good, 
the summaries should be all Passed with no Skipped or Failed tests.

Please build examples and .test scripts as you add features to the code.
  

Code documentation
------------------
is generated from embedded comments by Doxygen and starts at the 
file ./html/index.html.

After being unsatisfied with any of the available documentation
extraction systems, we've (i.e. PJ has) returned to Doxygen 
as the best available option.
The (lesser) evil of extra markup making code comments difficult
to read is offset by the ease with which one can browse the 
generated HTML documents.

The process of adding Doxygen markup to the comments 
has just been restarted so it's going to take some time 
for it to be comprehensive.
 
