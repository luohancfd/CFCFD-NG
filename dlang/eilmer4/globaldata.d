// globaldata.d
// Peter J. 2014-07-18 first cut.

module globaldata;

import block;
import flowstate;

// Each task/process look after a local "bag" of blocks 
// that may not be sequentially numbered.
Block[] allBlocks; // The array of Block objects, holding arrays of cells.
Block[] myBlocks;  // Local collection that we can iterate over.

// A place to store gas and flow properties for boundary conditions, etc.
static FlowState[] myFlowStates; 
