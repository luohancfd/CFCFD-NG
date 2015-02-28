// globaldata.d
// Peter J. 2014-07-18 first cut.

module globaldata;

import sblock;

// When we get around to implementing the MPI version of the code,
// each task/process look after a local "bag" of blocks 
// that may or may not be sequentially numbered.

static SBlock[] allBlocks; // The array of Block objects, holding arrays of cells.
static SBlock[] myBlocks;  // Local collection that we can iterate over.
