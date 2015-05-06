/**
 * globaldata.d
 *
 * Author: Peter J. and Rowan G.
 * Versions: 2014-07-18 : first cut.
 *           2015-22-04 : added containers for solid blocks
 */

module globaldata;

import globalconfig;
import sblock;
import ssolidblock;

// When we get around to implementing the MPI version of the code,
// each task/process look after a local "bag" of blocks 
// that may or may not be sequentially numbered.

static SBlock[] allBlocks; // The array of Block objects, holding arrays of cells.
static SBlock[] myBlocks;  // Local collection that we can iterate over.

static SSolidBlock[] allSolidBlocks;
static SSolidBlock[] mySolidBlocks;

// The current parallel code is based on having one SBlock per thread.
// We need to hava a dedicated set of configuration parameters for each thread so that
// there is no need to have memory barriers guarding their access.
static LocalConfig[] myConfig;
